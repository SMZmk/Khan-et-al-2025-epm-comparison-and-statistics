# Python script to generate SHAP saliency plots (Comparing UP vs DOWN Regulated Genes)

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import re
import h5py # Required for HDF5 files
import glob # For finding files with patterns

# --- Configuration ---
# *** Keep your original output directory or change as needed ***
output_dir_updown = "/home/ibg-4/Desktop/arab_env_2024i/studies/dCRE_results_AnemAsag20250130/results/shap/shap_plots_up_vs_down"

species_label_map = {
    'Anem': '$A.\ nemorensis$',
    'Asag': '$A.\ sagittata$'
}
condition_label_map = {
    'D0h': 'Control',
    'DWh': 'Wilting',
    # 'DSh': 'Recovery' # Add back if recovery data/models exist and are needed
}

# Define colors for UP and DOWN regulated lines
up_color = "#D55E00" # Vermillion/Orange
down_color = "#0072B2" # Blue

try:
    sns.set_theme(style="ticks", context="paper")
    plt.rcParams.update({
        'figure.dpi': 100, 'savefig.dpi': 300, 'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 11,
        'axes.labelsize': 11, 'axes.titlesize': 12, 'xtick.labelsize': 10,
        'ytick.labelsize': 9, 'legend.fontsize': 9, 'figure.titlesize': 14,
        'figure.constrained_layout.use': True, 'text.usetex': False
    })
except Exception as e:
    print(f"Error setting plotting parameters: {e}")

# --- File Path Configuration ---
shap_data_dir = "/home/ibg-4/Desktop/arab_env_2024i/studies/dCRE_results_AnemAsag20250130/results/shap"
dge_data_dir = "/home/ibg-4/Desktop/arab_env_2024i/runs/Arabis_nem_sag_drought_mRNA_counts"
window_size = 40
h5_file_pattern = "{spec}-{cond}-*_deepcre_interpret_*.h5"
meta_file_pattern = "{spec}-{cond}-*_deepcre_interpret_*_shap_meta.csv"
# Assumes DGE file compares Control vs Wilting (D0vW)
dge_file_pattern = "{spec}-D0vW-dge1_reg0.csv"


# --- Function to Load and Process SHAP Data (Unchanged) ---
def load_process_shap_data(species_code, condition_code, shap_dir):
    """Loads HDF5 SHAP scores and metadata for a given species and condition."""
    print(f"Processing SHAP data for: {species_code} - {condition_code}")
    h5_pattern = os.path.join(shap_dir, h5_file_pattern.format(spec=species_code, cond=condition_code))
    meta_pattern = os.path.join(shap_dir, meta_file_pattern.format(spec=species_code, cond=condition_code))
    h5_files = sorted(glob.glob(h5_pattern))
    meta_files = sorted(glob.glob(meta_pattern))
    if not h5_files or not meta_files:
        print(f"  Error: No H5 or Meta file found for {species_code}-{condition_code}.")
        return None, None
    h5_file_path = h5_files[0]
    h5_base = os.path.basename(h5_file_path).replace('.h5', '')
    match = re.match(r"(.*?_deepcre_interpret_\d+)_.*", h5_base)
    h5_core_id = match.group(1) if match else h5_base
    meta_file_path = None
    for mf in meta_files:
        if h5_core_id in mf: meta_file_path = mf; break
    if meta_file_path is None: meta_file_path = meta_files[0]; print(f"  Warning: Meta file match failed. Using {meta_file_path}.")
    print(f"  Using H5 file: {os.path.basename(h5_file_path)}")
    print(f"  Using Meta file: {os.path.basename(meta_file_path)}")

    try:
        with h5py.File(h5_file_path, 'r') as hf:
            if 'contrib_scores' not in hf: print(f"  Error: 'contrib_scores' not found."); return None, None
            contrib_scores = hf['contrib_scores'][:]
        df_meta = pd.read_csv(meta_file_path)
        if "gene_ids" in df_meta.columns and "gene_id" not in df_meta.columns: df_meta = df_meta.rename(columns={"gene_ids": "gene_id"})
        if "gene_id" not in df_meta.columns: print(f"  Error: No 'gene_id' column found."); return None, None

        if contrib_scores.ndim != 3 or contrib_scores.shape[2] != 4: print(f"  Error: Unexpected shape {contrib_scores.shape}. Expected (num_seq, seq_len, 4)."); return None, None
        num_seq_h5 = contrib_scores.shape[0]
        num_seq_meta = len(df_meta)
        if num_seq_h5 != num_seq_meta: print(f"  Warning: Mismatch H5 sequences ({num_seq_h5}) vs metadata ({num_seq_meta})."); return None, None # Changed to warning, allow partial processing if desired, but safer to return None

        # Aggregate scores (sum over nucleotides) - Keep original sign
        agg_scores_genes_pos = contrib_scores.sum(axis=2)

        print(f"  Successfully processed. Aggregated scores shape: {agg_scores_genes_pos.shape}")
        return agg_scores_genes_pos, df_meta

    except Exception as e:
        print(f"  Error reading or processing file {h5_file_path} or {meta_file_path}: {e}")
        return None, None

# --- Function to Calculate Positive and Negative Means (Unchanged) ---
def calculate_pos_neg_means(scores_matrix):
    """Calculates the mean of positive and negative scores separately per position."""
    if scores_matrix is None or scores_matrix.shape[0] == 0:
        return None, None

    num_genes, num_positions = scores_matrix.shape
    mean_pos_scores = np.zeros(num_positions)
    mean_neg_scores = np.zeros(num_positions)

    # Calculate mean where scores > 0
    pos_scores = np.where(scores_matrix > 0, scores_matrix, np.nan)
    pos_counts = np.sum(~np.isnan(pos_scores), axis=0)
    pos_sums = np.nansum(pos_scores, axis=0)
    mean_pos_scores = np.divide(pos_sums, pos_counts, where=pos_counts > 0, out=np.zeros_like(pos_sums))

    # Calculate mean where scores < 0
    neg_scores = np.where(scores_matrix < 0, scores_matrix, np.nan)
    neg_counts = np.sum(~np.isnan(neg_scores), axis=0)
    neg_sums = np.nansum(neg_scores, axis=0)
    mean_neg_scores = np.divide(neg_sums, neg_counts, where=neg_counts > 0, out=np.zeros_like(neg_sums))

    return mean_pos_scores, mean_neg_scores


# --- Function to get rolling mean scores (Pos and Neg), optionally filtered (Unchanged) ---
def get_pos_neg_rolling_means(agg_scores, meta_df, gene_ids_to_keep=None):
    """Calculates pos/neg means and rolling means, optionally filtering genes."""
    if agg_scores is None or meta_df is None:
        print("  Warning: Missing scores or metadata for rolling mean calculation.")
        return None, None

    filtered_scores = None
    if gene_ids_to_keep is not None:
        # Ensure gene_ids_to_keep is a set for faster lookup
        gene_ids_set = set(gene_ids_to_keep)
        if not gene_ids_set:
             print("  Warning: Provided gene list for filtering is empty.")
             return None, None

        meta_df_reset = meta_df.reset_index() # Ensure proper indexing if meta_df had non-standard index
        # Use boolean indexing which is generally safer and clearer
        keep_mask = meta_df_reset['gene_id'].isin(gene_ids_set)
        indices = meta_df_reset.index[keep_mask].tolist()

        if len(indices) == 0:
            print(f"  Warning: No matching genes found in metadata for the {len(gene_ids_set)} provided gene IDs.")
            return None, None

        # Check bounds before indexing agg_scores
        valid_indices = [idx for idx in indices if idx < agg_scores.shape[0]]
        if len(valid_indices) != len(indices):
             print(f"  Warning: Some indices ({len(indices) - len(valid_indices)}) were out of bounds for the scores matrix shape {agg_scores.shape}. Using {len(valid_indices)} valid indices.")
             if not valid_indices:
                 return None, None # No valid indices remain
             indices = valid_indices

        filtered_scores = agg_scores[indices, :]
        if filtered_scores.shape[0] == 0:
             print("  Warning: Filtering resulted in zero genes after index validation.")
             return None, None
        print(f"  Filtered scores shape for rolling mean: {filtered_scores.shape} (using {len(gene_ids_set)} target IDs)")
        scores_to_process = filtered_scores
    else:
        # Use all scores
        scores_to_process = agg_scores
        print(f"  Using all scores for rolling mean. Shape: {scores_to_process.shape}")


    # Calculate separate means
    mean_pos, mean_neg = calculate_pos_neg_means(scores_to_process)
    if mean_pos is None or mean_neg is None:
        print("  Warning: Calculation of pos/neg means failed.")
        return None, None

    # Apply rolling mean to each
    rolling_mean_pos = pd.Series(mean_pos).rolling(window=window_size, center=True, min_periods=1).mean()
    rolling_mean_neg = pd.Series(mean_neg).rolling(window=window_size, center=True, min_periods=1).mean()

    return rolling_mean_pos, rolling_mean_neg


# --- MODIFIED Function to Plot UP vs DOWN Saliency Comparison ---
def plot_up_down_saliency_comparison(
    series_up_pos, series_up_neg, series_down_pos, series_down_neg,
    color_up, color_down, title, filename, output_dir):
    """Plots positive and negative rolling mean series comparing UP vs DOWN regulated genes."""
    fig, ax = plt.subplots(figsize=(10, 5)) # Slightly taller figure

    # Check if all data is present for *both* UP and DOWN sets
    if series_up_pos is None or series_up_neg is None:
        print(f"Skipping plot '{title}': Missing data for UP regulated genes.")
        plt.close(fig)
        return
    if series_down_pos is None or series_down_neg is None:
        print(f"Skipping plot '{title}': Missing data for DOWN regulated genes.")
        plt.close(fig)
        return

    # Determine positions based on one of the series (assuming they all have the same length)
    positions = np.arange(len(series_up_pos)) - (len(series_up_pos) // 2) # Center positions around 0

    # Plotting - Use different linestyles for pos/neg
    label_up = "UP Regulated"
    label_down = "DOWN Regulated"

    # UP Regulated Genes
    ax.plot(positions, series_up_pos.values, label=f"{label_up} (Positive Mean)", color=color_up, linewidth=1.5, linestyle='-')
    ax.plot(positions, series_up_neg.values, label=f"{label_up} (Negative Mean)", color=color_up, linewidth=1.5, linestyle=':') # Dotted for neg

    # DOWN Regulated Genes
    ax.plot(positions, series_down_pos.values, label=f"{label_down} (Positive Mean)", color=color_down, linewidth=1.5, linestyle='-')
    ax.plot(positions, series_down_neg.values, label=f"{label_down} (Negative Mean)", color=color_down, linewidth=1.5, linestyle=':') # Dotted for neg

    ax.axhline(0, color='grey', linestyle='-', linewidth=0.8) # Solid line at 0
    ax.set_ylabel("Mean SHAP Score (Rolling Mean)")
    ax.set_title(title)
    ax.legend(fontsize=plt.rcParams['legend.fontsize'] * 0.9)

    ax.axvline(0, color='dimgrey', linestyle='--', linewidth=0.8) # Line at feature center
    ax.set_xlabel("Position relative to feature center (bp)")

    os.makedirs(output_dir, exist_ok=True)
    full_filename_pdf = os.path.join(output_dir, f"{filename}.pdf")
    full_filename_png = os.path.join(output_dir, f"{filename}_600dpi.png")
    try:
        fig.savefig(full_filename_pdf, bbox_inches='tight', pad_inches=0.1)
        fig.savefig(full_filename_png, dpi=600, bbox_inches='tight', pad_inches=0.1)
        print(f"  Saved plot: {os.path.basename(full_filename_pdf)}")
        # print(f"  Saved plot: {os.path.basename(full_filename_png)}") # Optional: uncomment if you want both messages
    except Exception as e:
        print(f"  Error saving plot {filename}: {e}")

    plt.close(fig)


# --- MODIFIED Main Execution Logic ---
def main():
    os.makedirs(output_dir_updown, exist_ok=True)
    print(f"Created output directory for UP vs DOWN plots: {output_dir_updown}")

    # --- Load SHAP Data (Raw Scores and Meta) ---
    raw_data_cache = {}
    species_to_process = ['Anem', 'Asag']
    conditions_to_process = ['D0h', 'DWh'] # Add 'DSh' if needed

    print("\n--- Loading and Processing SHAP Data ---")
    for spec in species_to_process:
        for cond in conditions_to_process:
            key = f"{spec}_{cond}"
            agg_scores, meta_df = load_process_shap_data(spec, cond, shap_data_dir)
            if agg_scores is not None and meta_df is not None:
                raw_data_cache[key] = {'scores': agg_scores, 'meta': meta_df}
            else:
                print(f"Could not load or process data for {key}.")
                raw_data_cache[key] = None # Store None to indicate failure

    # --- Load DGE Data ---
    print("\n--- Loading DGE Data ---")
    dge_data = {}
    for spec in species_to_process:
        dge_file = os.path.join(dge_data_dir, dge_file_pattern.format(spec=spec))
        try:
            dge_df = pd.read_csv(dge_file)
            # Basic check for required columns
            if 'gene_id' in dge_df.columns and 'regulation' in dge_df.columns:
                 # Remove potential duplicates just in case
                dge_df = dge_df.drop_duplicates(subset=['gene_id'])
                dge_data[spec] = dge_df
                print(f"  Loaded DGE data for {spec} from {os.path.basename(dge_file)}")
            else:
                print(f"  Error: Missing 'gene_id' or 'regulation' columns in {dge_file}");
                dge_data[spec] = None
        except FileNotFoundError:
            print(f"  Error: DGE file not found: {dge_file}");
            dge_data[spec] = None
        except Exception as e:
            print(f"  Error reading DGE file {dge_file}: {e}");
            dge_data[spec] = None

    # --- Extract UP/DO Gene Lists ---
    gene_lists = {'Anem': {'UP': [], 'DO': []}, 'Asag': {'UP': [], 'DO': []}}
    print("\n--- Extracting Gene Lists ---")
    for spec, df in dge_data.items():
        if df is not None:
            # Ensure regulation column has expected values
            up_genes = df[df['regulation'] == 'UP']['gene_id'].unique().tolist()
            down_genes = df[df['regulation'] == 'DO']['gene_id'].unique().tolist()
            gene_lists[spec]['UP'] = up_genes
            gene_lists[spec]['DO'] = down_genes
            print(f"  {spec}: Found {len(up_genes)} unique UP genes, {len(down_genes)} unique DO genes (from D0vW).")
        else:
             print(f"  {spec}: No DGE data loaded, cannot extract gene lists.")


    # --- Generate Plots: Iterate through Species and Conditions ---
    print("\n--- Generating UP vs DOWN Saliency Plots ---")
    for spec in species_to_process:
        spec_label = species_label_map[spec]
        up_gene_list = gene_lists[spec].get('UP')
        down_gene_list = gene_lists[spec].get('DO')

        # Check if we have both UP and DOWN lists for this species
        if not up_gene_list:
            print(f"Skipping plots for {spec_label}: Missing UP regulated gene list.")
            continue
        if not down_gene_list:
            print(f"Skipping plots for {spec_label}: Missing DOWN regulated gene list.")
            continue
        if not up_gene_list and not down_gene_list:
             print(f"Skipping plots for {spec_label}: Missing both UP and DOWN regulated gene lists.")
             continue


        for cond_code in conditions_to_process:
            cond_label = condition_label_map[cond_code]
            data_key = f"{spec}_{cond_code}"
            shap_data = raw_data_cache.get(data_key) # Retrieve cached data

            # Check if SHAP data was loaded successfully for this combo
            if not shap_data:
                print(f"Skipping plot for {spec_label} - {cond_label}: SHAP data not loaded.")
                continue

            print(f"\n--- Generating Plot: {spec_label} - {cond_label} (UP vs DOWN) ---")

            current_scores = shap_data.get('scores')
            current_meta = shap_data.get('meta')

            # Get rolling means for UP genes using the current condition's model
            print("  Calculating rolling means for UP genes...")
            up_pos, up_neg = get_pos_neg_rolling_means(
                current_scores,
                current_meta,
                gene_ids_to_keep=up_gene_list
            )

            # Get rolling means for DOWN genes using the *same* current condition's model
            print("  Calculating rolling means for DOWN genes...")
            down_pos, down_neg = get_pos_neg_rolling_means(
                current_scores,
                current_meta,
                gene_ids_to_keep=down_gene_list
            )

            # Plot the UP vs DOWN comparison for this species and condition
            plot_up_down_saliency_comparison(
                series_up_pos=up_pos, series_up_neg=up_neg,
                series_down_pos=down_pos, series_down_neg=down_neg,
                color_up=up_color, color_down=down_color,
                title=f"Mean Saliency ({spec_label}, {cond_label} Model) - UP vs. DOWN Genes",
                filename=f"saliency_posneg_UPvsDOWN_{spec}_{cond_code}",
                output_dir=output_dir_updown # Use the new output directory
            )

    print(f"\nProcessing finished. UP vs DOWN plots saved in {output_dir_updown}")

if __name__ == "__main__":
    main()
