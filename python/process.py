#!/usr/bin/env python3

"""
Variant annotation comparison analysis

This script:

1. Reads all sheets from an Excel workbook
   (one chromosome per sheet)

2. Calculates:
   - global concordance
   - per-category concordance
   - reclassification rate
   - clinically relevant upgrades
   - clinically relevant downgrades

3. Generates:
   - raw heatmaps
   - normalized heatmaps
   - stacked barplots
   - chromosome summary heatmap
   - clinical upgrade plots

IMPORTANT:
The category "other" is EXCLUDED from ALL plots,
but remains included in the statistical calculations.

Rows = OLD annotation
Columns = NEW annotation
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

INPUT_EXCEL = "CHROMS_DATOS.xlsx"
OUTPUT_DIR = "variant_analysis_output"
os.makedirs(OUTPUT_DIR, exist_ok=True)


# CATEGORY ORDER
categories = [
    "unclassified",
    "other",
    "benign",
    "likely benign",
    "vus",
    "conflicting",
    "likely pathogenic",
    "pathogenic"
]

plot_categories = [
    "unclassified",
    "benign",
    "likely benign",
    "vus",
    "conflicting",
    "likely pathogenic",
    "pathogenic"
]

clinical_categories = [
    "benign",
    "likely benign",
    "vus",
    "conflicting",
    "likely pathogenic",
    "pathogenic"
]


xls = pd.ExcelFile(INPUT_EXCEL)
summary_rows = []


def calculate_concordance(df):
    total = df.values.sum()
    diagonal = np.trace(df.values)
    if total == 0:
        return np.nan
    return diagonal / total

def calculate_row_concordance(df):
    row_conc = {}
    for cat in df.index:
        row_total = df.loc[cat].sum()
        if row_total == 0:
            row_conc[cat] = np.nan
        else:
            row_conc[cat] = (
                df.loc[cat, cat] / row_total
            )
    return row_conc

def calculate_clinical_upgrades(df):
    upgrades = {}
    upgrades["to_likely_pathogenic"] = (
        df["likely pathogenic"].sum()
        - df.loc[
            "likely pathogenic",
            "likely pathogenic"
        ]
    )

    upgrades["to_pathogenic"] = (
        df["pathogenic"].sum()
        - df.loc[
            "pathogenic",
            "pathogenic"
        ]
    )

    upgrades["vus_to_likely_pathogenic"] = (
        df.loc["vus", "likely pathogenic"]
    )

    upgrades["vus_to_pathogenic"] = (
        df.loc["vus", "pathogenic"]
    )

    upgrades["benign_to_pathogenic"] = (
        df.loc["benign", "pathogenic"]
    )

    upgrades["benign_to_likely_pathogenic"] = (
        df.loc["benign", "likely pathogenic"]
    )

    upgrades["likely_benign_to_pathogenic"] = (
        df.loc["likely benign", "pathogenic"]
    )

    return upgrades


def calculate_clinical_downgrades(df):
    downgrades = {}
    downgrades["pathogenic_to_vus"] = (
        df.loc["pathogenic", "vus"]
    )
    downgrades["pathogenic_to_benign"] = (
        df.loc["pathogenic", "benign"]
    )
    downgrades["likely_pathogenic_to_vus"] = (
        df.loc["likely pathogenic", "vus"]
    )
    downgrades["likely_pathogenic_to_benign"] = (
        df.loc["likely pathogenic", "benign"]
    )
    return downgrades


for sheet in xls.sheet_names:
    print(f"Processing {sheet}")
    df = pd.read_excel(
        xls,
        sheet_name=sheet,
        index_col=0
    )

    df.index = (
        df.index
        .astype(str)
        .str.strip()
        .str.lower()
    )

    df.columns = (
        df.columns
        .astype(str)
        .str.strip()
        .str.lower()
    )

    categories_clean = [
        c.strip().lower()
        for c in categories
    ]

    df = df.reindex(
        index=categories_clean,
        columns=categories_clean
    ).fillna(0)

    df = (
        df.apply(
            pd.to_numeric,
            errors="coerce"
        )
        .fillna(0)
        .astype(int)
    )

    if df.values.sum() == 0:
        print(f"Skipping empty chromosome: {sheet}")
        continue

    total_variants = df.values.sum()
    concordance = calculate_concordance(df)
    reclassification_rate = 1 - concordance
    row_concordance = calculate_row_concordance(df)
    upgrades = calculate_clinical_upgrades(df)
    downgrades = calculate_clinical_downgrades(df)

    row = {
        "chromosome": sheet,
        "total_variants": total_variants,
        "concordance": concordance,
        "reclassification_rate": reclassification_rate
    }

    for k, v in row_concordance.items():
        row[f"{k}_stability"] = v

    row.update(upgrades)
    row.update(downgrades)
    summary_rows.append(row)
    plot_df = df.loc[
        plot_categories,
        plot_categories
    ]

    ##REPRESENTACIONES GRÑAFICAS
    #HEATMAP
    color_df = plot_df.copy().astype(float)
    color_df[color_df == 0] = 0.5
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        color_df,
        cmap="viridis",
        norm=LogNorm(
            vmin=1,
            vmax=color_df.values.max()
        ),
        annot=plot_df,      # show original values
        fmt=".0f",
        linewidths=0.5,
        linecolor="white"
    )

    plt.title(f"Transiciones ACMG en el {sheet}")
    plt.xlabel("Última versión de la anotación")
    plt.ylabel("Versión anterior de la anotación")
    plt.tight_layout()
    plt.savefig(
        os.path.join(
            OUTPUT_DIR,
            f"{sheet}_raw_heatmap.png"
        ),
        dpi=300
    )

    plt.close()

    #HEATMAP NORMALIZADO
    row_sums = plot_df.sum(axis=1)
    norm_df = plot_df.div(
        row_sums.replace(0, np.nan),
        axis=0
    )

    if not np.all(np.isnan(norm_df.values)):
        plt.figure(figsize=(10, 8))
        sns.heatmap(
            norm_df,
            cmap="magma",
            annot=True,
            fmt=".2%"
        )
        plt.title(
            f"Transiciones ACMG en el {sheet} normalizadas"
        )
        plt.xlabel("Última versión de la anotación")
        plt.ylabel("Versión anterior de la anotación")
        plt.tight_layout()
        plt.savefig(
            os.path.join(
                OUTPUT_DIR,
                f"{sheet}_normalized_heatmap.png"
            ),
            dpi=300
        )
        plt.close()

    #BARPLOT
    plt.figure(figsize=(12, 8))
    norm_df.plot(
        kind="bar",
        stacked=True
    )
    plt.title(
        f"Transiciones ACMG en el {sheet}"
    )
    plt.ylabel("Porcentaje")
    plt.xlabel("Versión anterior de la anotación")
    plt.legend(
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        title= "Última versión de la anotación"
    )
    plt.tight_layout()
    plt.savefig(
        os.path.join(
            OUTPUT_DIR,
            f"{sheet}_stacked_barplot.png"
        ),
        dpi=300
    )
    plt.close()

    #HEATMAP CLASIFICADAS
    clinical_df = plot_df.loc[
        clinical_categories,
        clinical_categories
    ]
    clinical_norm = clinical_df.div(
        clinical_df.sum(axis=1).replace(0, np.nan),
        axis=0
    )

    if not np.all(np.isnan(clinical_norm.values)):
        plt.figure(figsize=(8, 6))
        sns.heatmap(
            clinical_norm,
            cmap="coolwarm",
            annot=True,
            fmt=".2%"
        )
        plt.title(f"Transiciones ACMG en el {sheet} para variantes clasificadas"
        )
        plt.xlabel("Última versión de la anotación")
        plt.ylabel("Versión anterior de la anotación")
        plt.tight_layout()
        plt.savefig(
            os.path.join(
                OUTPUT_DIR,
                f"{sheet}_clinical_heatmap.png"
            ),
            dpi=300
        )
        plt.close()


#Tabla resumen
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "summary_statistics.csv"
    ),
    index=False
)
summary_df.to_excel(
    os.path.join(
        OUTPUT_DIR,
        "summary_statistics.xlsx"
    ),
    index=False
)

#Heatmap resumen
summary_heatmap = summary_df[
    [
        "concordance",
        "reclassification_rate",
        "vus_stability",
        "pathogenic_stability",
        "to_pathogenic",
        "to_likely_pathogenic"
    ]
].copy()
summary_heatmap.index = (
    summary_df["chromosome"]
)
plt.figure(figsize=(12, 10))
sns.heatmap(
    summary_heatmap,
    annot=True,
    cmap="RdYlBu_r"
)
plt.title(
    "Resumen de las métricas analizadas por cromosoma"
)
plt.tight_layout()
plt.savefig(
    os.path.join(
        OUTPUT_DIR,
        "chromosome_summary_heatmap.png"
    ),
    dpi=300
)
plt.close()


#Barplot variantes patogénicas
upgrade_cols = [
    "to_likely_pathogenic",
    "to_pathogenic"
]

upgrade_df = summary_df[
    ["chromosome"] + upgrade_cols
].set_index("chromosome")

upgrade_df.plot(
    kind="bar",
    figsize=(14, 8)
)
plt.title(
    "Variantes que pasan a categorías patogénicas detectadas por cromosoma"
)
plt.ylabel("Número de variantes")
plt.xlabel("Cromosoma")
plt.tight_layout()
plt.savefig(
    os.path.join(
        OUTPUT_DIR,
        "clinical_upgrades.png"
    ),
    dpi=300
)
plt.close()

print(
    f"Results saved in: {OUTPUT_DIR}"
)