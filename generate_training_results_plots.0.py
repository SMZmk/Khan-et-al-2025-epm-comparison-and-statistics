# Python script to generate plots (Version 7 - Final Labels)

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import re

# --- Configuration ---
data_dir = "/home/ibg-4/Desktop/arab_env_2024i/studies/results_arabis_2025d"
# Store plots in a final subdirectory
output_dir = os.path.join(data_dir, "plots_panel_1_v7_final_labels")
metrics_to_plot = ['accuracy', 'auROC', 'auPR']
file_pattern = re.compile(r"^(?P<Species>[A-Za-z]+)-(?P<Condition>D[0WS]h)-.*_train_models_(?P<Model_type>ssc|ssr)_.*\.csv$")

# --- Label Mappings ---
species_label_map = {
    'Anem': '$A.\ nemorensis$',
    'Asag': '$A.\ sagittata$'
}
condition_label_map = {
    'D0h': 'Control',
    'DWh': 'Wilting',
    'DSh': 'Recovery'
}

# --- Color Palette ---
plot_palette = "YlOrBr" # Earth tones

# --- Basic Plotting Configuration (Larger Fonts) ---
sns.set_theme(style="ticks", context="paper")
plt.rcParams.update({
    'figure.dpi': 100,
    'savefig.dpi': 300,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 11,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.titlesize': 14,
    'figure.constrained_layout.use': True,
    'text.usetex': False
})

# --- Data Loading and Processing ---
all_data = []
print(f"Looking for files in: {data_dir}")
try:
    filenames = os.listdir(data_dir)
except FileNotFoundError:
    print(f"Error: Data directory not found at {data_dir}")
    filenames = []

found_files = 0
for filename in filenames:
    match = file_pattern.match(filename)
    if match:
        try:
            # Keep original Condition code if needed, but map for plotting
            species = match.group('Species')
            condition_code = match.group('Condition') # Keep original code
            model_type = match.group('Model_type')
            file_path = os.path.join(data_dir, filename)
            df = pd.read_csv(file_path)
            if not all(col in df.columns for col in ['test'] + metrics_to_plot):
                 print(f"Warning: File {filename} missing required columns. Skipping.")
                 continue
            df['Species'] = species
            df['Condition'] = condition_code # Store original code if needed elsewhere
            # Add the new descriptive label column for plotting
            df['Treatment'] = condition_label_map.get(condition_code, condition_code)
            df['Model_type'] = model_type
            df['source_file'] = filename
            all_data.append(df)
            found_files += 1
        except Exception as e:
            print(f"Error processing file {filename}: {e}")
print(f"Found and processed {found_files} matching files.")


if not all_data:
    print("No data loaded. Cannot create plots.")
else:
    combined_df = pd.concat(all_data, ignore_index=True)
    for metric in metrics_to_plot:
        combined_df[metric] = pd.to_numeric(combined_df[metric], errors='coerce')
    # Drop rows where metric conversion failed OR mapping failed (unlikely with .get)
    combined_df = combined_df.dropna(subset=metrics_to_plot + ['Treatment'])


    print("\nCombined DataFrame head (with Treatment label):")
    print(combined_df.head())

    # --- Create Output Directory ---
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nCreated output directory: {output_dir}")

    # --- Generate Plots ---
    species_order = sorted(combined_df['Species'].unique())
    # Use the new descriptive labels for ordering
    treatment_order = ['Control', 'Wilting', 'Recovery']
    model_order = ['ssc', 'ssr']

    # Recalculate n_counts using the NEW 'Treatment' column for consistency in plot annotations
    n_counts = combined_df.groupby(['Species', 'Treatment', 'Model_type']).size().reset_index(name='n')
    print("\nCounts (n) per group (using Treatment labels):")
    print(n_counts)


    for metric in metrics_to_plot:
        print(f"\nGenerating plot for: {metric}")

        # 1. Create the main bar plot using 'Treatment' for hue
        g = sns.catplot(
            data=combined_df,
            kind='bar',
            x='Species',
            y=metric,
            hue='Treatment',      # <<< Use Treatment column for hue
            hue_order=treatment_order, # <<< Use new order
            col='Model_type',
            order=species_order,
            col_order=model_order,
            errorbar='sd',
            capsize=0.1,
            palette=plot_palette,
            legend_out=True,
            sharey=True,
        )

        # 2. Iterate through axes for overlays and labels
        for ax_idx, ax in enumerate(g.axes.flat):
            current_model_type = model_order[ax_idx % len(model_order)]
            # Filter data for the current facet's model type
            # No need to filter by Treatment here, stripplot handles hue
            facet_data = combined_df[combined_df['Model_type'] == current_model_type]

            # 3. Overlay individual data points using 'Treatment' for hue
            sns.stripplot(
                data=facet_data, x='Species', y=metric,
                hue='Treatment',      # <<< Use Treatment column for hue
                hue_order=treatment_order, # <<< Use new order
                order=species_order,
                marker='x', size=4, color='.2', alpha=0.7,
                jitter=True, dodge=True, legend=False, ax=ax
            )

            # --- Add 'n' Annotations (matching using Treatment labels) ---
            num_hues = len(treatment_order)
            for i, bar in enumerate(ax.patches):
                 loc_idx = i // num_hues
                 hue_idx = i % num_hues
                 if loc_idx < len(species_order):
                    current_species = species_order[loc_idx]
                    current_treatment = treatment_order[hue_idx] # Get label from order

                    # Find 'n' using the recalculated n_counts DataFrame
                    n_val_row = n_counts[
                        (n_counts['Species'] == current_species) &
                        (n_counts['Treatment'] == current_treatment) & # Match using Treatment label
                        (n_counts['Model_type'] == current_model_type)
                    ]

                    if not n_val_row.empty:
                        n_val = n_val_row['n'].iloc[0]
                        bar_height = bar.get_height()
                        # Use original combined_df for std dev calculation, filtering by treatment label
                        group_data = combined_df[
                            (combined_df['Species'] == current_species) &
                            (combined_df['Treatment'] == current_treatment) & # Filter using label
                            (combined_df['Model_type'] == current_model_type)
                        ]
                        std_dev = group_data[metric].std() if not group_data.empty else 0

                        error_bar_top = bar_height + std_dev
                        ylim_low, ylim_high = (0.5, 1.0)
                        text_y = error_bar_top + (ylim_high - ylim_low) * 0.03
                        text_y = min(text_y, ylim_high * 0.98)
                        text_y = max(text_y, ylim_low + (ylim_high - ylim_low) * 0.02)
                        ax.text(bar.get_x() + bar.get_width() / 2., text_y, f'n={n_val}',
                                ha='center', va='bottom', fontsize=7, color='black')


            # 4. Update X-tick labels with italics
            current_labels = [item.get_text() for item in ax.get_xticklabels()]
            new_labels = [species_label_map.get(label, label) for label in current_labels]
            ax.set_xticklabels(new_labels)
            ax.set_xlabel("")

        # Set Y-axis Limits
        print(f"  Setting Y-axis limits to (0.5, 1.0) for {metric}")
        g.set(ylim=(0.5, 1.0))

        # Improve overall plot aesthetics
        g.set_axis_labels("Species", f"Mean {metric.upper()}")
        g.set_titles("Model: {col_name}")
        g.fig.suptitle(f'Model Performance Comparison: {metric.upper()}', y=1.03)
        try:
           # Update legend title to reflect the new labels
           g.legend.set_title("Treatment") # <<< Changed legend title
           legend_title_obj = g.legend.get_title()
           legend_title_obj.set_fontsize(plt.rcParams['legend.fontsize'])
           sns.move_legend(g, "upper left", bbox_to_anchor=(1, 0.9))
        except AttributeError:
           print("Could not modify legend attributes directly.")

        try:
            # Rely on constrained_layout
            pass
        except ValueError:
             print("Warning: Layout adjustment failed.")

        # --- Save the plot ---
        plot_filename_base = f"{metric}_comparison_final" # Simplified filename
        pdf_path = os.path.join(output_dir, f"{plot_filename_base}.pdf")
        png_path = os.path.join(output_dir, f"{plot_filename_base}_600dpi.png")

        try:
            g.fig.savefig(pdf_path, bbox_inches='tight', pad_inches=0.1)
            print(f"  Saved PDF: {pdf_path}")
            g.fig.savefig(png_path, dpi=600, bbox_inches='tight', pad_inches=0.1)
            print(f"  Saved PNG (600 DPI): {png_path}")
        except Exception as e:
            print(f"  Error saving plot for {metric}: {e}")

        plt.close(g.fig)

    print("\nFinished generating all plots.")
