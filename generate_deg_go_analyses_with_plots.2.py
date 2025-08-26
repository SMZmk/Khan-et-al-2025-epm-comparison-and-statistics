# Import necessary libraries
import pandas as pd
import re
import os
import sys
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from collections import defaultdict
import math # For -log10

# --- Configuration ---
start_df_path = "merged_ortho_quadrants_go_enriched.tsv"
motif_dir = "/home/ibg-4/Desktop/arab_env_2024i/studies/AnemAsag_EPMclu"
cluster_dir = "/home/ibg-4/Desktop/arab_env_2024i/studies/AnemAsag_EPM_set_analyses_K5K9"
anem_motif_filename = "2025bAnemSWAnem_AnemE-04_SW2025b_gene_none-q1q9.csv"
asag_motif_filename = "2025dAsagSWAsagSW_E-04_2025d_gene_none-q1q9.csv"
cluster_filename = "EPM_cluster_assignments_k25.csv"
anem_motif_path = os.path.join(motif_dir, anem_motif_filename)
asag_motif_path = os.path.join(motif_dir, asag_motif_filename)
cluster_path = os.path.join(cluster_dir, cluster_filename)
output_dir = "quadrant_go_motif_networks_stricter" # New output dir name
os.makedirs(output_dir, exist_ok=True)

# --- Analysis Parameters ---
MIN_GENES_FOR_NODE = 5 # Stricter: Minimum A.thaliana genes for a GO term node
GO_PVAL_THRESHOLD = 0.01 # Stricter: Max corrected p-value for GO term node

# --- Helper Functions ---
def remove_version(gene_id):
    if pd.isna(gene_id): return None
    return re.sub(r'\.\d+$', '', str(gene_id))

def parse_go_term_string_with_pval(go_string):
    """Parses 'GO:ID (Name) [p_corr=value]' -> (GO_ID, Name, p_value)."""
    if pd.isna(go_string) or not isinstance(go_string, str) or not go_string.strip():
        return None, None, None
    # Match GO ID, Name, and p_corr value
    match = re.match(r'(GO:\d+)\s+\((.*?)\)(?:\s+\[p_corr=([\d.eE\-+]+)\])?', go_string)
    if match:
        go_id = match.group(1)
        go_name = match.group(2)
        p_val_str = match.group(3)
        p_value = None
        if p_val_str:
            try:
                p_value = float(p_val_str)
            except ValueError:
                p_value = None # Handle if p-value format is unexpected
        return go_id, go_name, p_value
    else:
        # Fallback if only GO ID is present
        match_simple = re.match(r'(GO:\d+)', go_string)
        if match_simple:
            return match_simple.group(1), match_simple.group(1), None # No p-value info
        else:
            return None, None, None

# --- 1. Load Start DataFrame ---
print(f"Loading base data from {start_df_path}...")
try:
    start_DF = pd.read_csv(start_df_path, sep='\t')
    print(f"Loaded {len(start_DF)} rows.")
except FileNotFoundError: print(f"Error: File not found {start_df_path}."); sys.exit(1)
except Exception as e: print(f"Error loading {start_df_path}: {e}"); sys.exit(1)

# --- 2. Load/Process Motif Data ---
# (Same as previous version - load, clean, filter)
# ... [omitted for brevity] ...
print(f"Loading/Processing Anem motif data from {anem_motif_path}...")
anem_motif_df = None
try:
    anem_motif_df = pd.read_csv(anem_motif_path)[['loc_ID', 'motif']].rename(columns={'loc_ID': 'gene_id'})
    anem_motif_df['gene_id_cleaned'] = anem_motif_df['gene_id'].apply(remove_version)
    anem_motif_df = anem_motif_df.dropna(subset=['gene_id_cleaned', 'motif'])
    print(f"  Loaded {len(anem_motif_df)} Anem associations.")
except Exception as e: print(f"  Error Anem motifs: {e}")
print(f"Loading/Processing Asag motif data from {asag_motif_path}...")
asag_motif_df = None
try:
    asag_motif_df = pd.read_csv(asag_motif_path)[['loc_ID', 'motif']].rename(columns={'loc_ID': 'gene_id'})
    asag_motif_df['gene_id_cleaned'] = asag_motif_df['gene_id'].apply(remove_version)
    asag_motif_df = asag_motif_df.dropna(subset=['gene_id_cleaned', 'motif'])
    print(f"  Loaded {len(asag_motif_df)} Asag associations.")
except Exception as e: print(f"  Error Asag motifs: {e}")
anem_ids_in_start_df = set(start_DF['Anemorosa_cleaned_ID'].dropna())
asag_ids_in_start_df = set(start_DF['Asagittata_cleaned_ID'].dropna())
if anem_motif_df is not None:
    anem_motif_df = anem_motif_df[anem_motif_df['gene_id_cleaned'].isin(anem_ids_in_start_df)]
    print(f"  Filtered Anem motifs: {len(anem_motif_df)} rows.")
if asag_motif_df is not None:
    asag_motif_df = asag_motif_df[asag_motif_df['gene_id_cleaned'].isin(asag_ids_in_start_df)]
    print(f"  Filtered Asag motifs: {len(asag_motif_df)} rows.")

# --- 3. Load/Process Cluster Data ---
# (Same as previous version - load, trim epm, create mapping)
# ... [omitted for brevity] ...
print(f"Loading/Processing cluster data from {cluster_path}...")
cluster_df = None
try:
    cluster_df = pd.read_csv(cluster_path)[['epm', 'cluster']]
    cluster_df['epm_trimmed'] = cluster_df['epm'].astype(str).str[:17]
    cluster_df = cluster_df[['epm_trimmed', 'cluster']].drop_duplicates(subset=['epm_trimmed'], keep='first')
    print(f"  Loaded {len(cluster_df)} unique trimmed epm-cluster assignments.")
except Exception as e: print(f"  Error Cluster data: {e}")

# --- 4. Add Cluster Info & Create Mappings ---
# (Same as previous version)
# ... [omitted for brevity] ...
anem_gene_to_motifclust = defaultdict(set); asag_gene_to_motifclust = defaultdict(set)
if anem_motif_df is not None and cluster_df is not None:
    anem_motifs_clust = pd.merge(anem_motif_df, cluster_df, left_on='motif', right_on='epm_trimmed', how='left')
    anem_motifs_clust['cluster'] = anem_motifs_clust['cluster'].fillna(0).astype(int)
    grouped = anem_motifs_clust.groupby('gene_id_cleaned')[['motif', 'cluster']]
    for name, group in grouped: anem_gene_to_motifclust[name] = set(map(tuple, group.drop_duplicates().values))
    print(f"  Created Anem gene -> (motif, cluster) map for {len(anem_gene_to_motifclust)} genes.")
if asag_motif_df is not None and cluster_df is not None:
    asag_motifs_clust = pd.merge(asag_motif_df, cluster_df, left_on='motif', right_on='epm_trimmed', how='left')
    asag_motifs_clust['cluster'] = asag_motifs_clust['cluster'].fillna(0).astype(int)
    grouped = asag_motifs_clust.groupby('gene_id_cleaned')[['motif', 'cluster']]
    for name, group in grouped: asag_gene_to_motifclust[name] = set(map(tuple, group.drop_duplicates().values))
    print(f"  Created Asag gene -> (motif, cluster) map for {len(asag_gene_to_motifclust)} genes.")


# --- 5. Prepare Mappings and GO Info from start_DF ---
print("Processing GO/Gene/Ortholog info from start_DF...")
ath_to_anem = start_DF.dropna(subset=['Athaliana_cleaned_ID', 'Anemorosa_cleaned_ID']).groupby('Athaliana_cleaned_ID')['Anemorosa_cleaned_ID'].apply(set)
ath_to_asag = start_DF.dropna(subset=['Athaliana_cleaned_ID', 'Asagittata_cleaned_ID']).groupby('Athaliana_cleaned_ID')['Asagittata_cleaned_ID'].apply(set)
go_info_per_quadrant = defaultdict(list) # Store list of (GO_ID, GO_Name, P_value) tuples
go_term_to_genes_per_quad = defaultdict(lambda: defaultdict(set)) # Stores GO_ID -> {gene set} per quadrant

if 'Quadrant_Enriched_GOs' not in start_DF.columns: start_DF['Quadrant_Enriched_GOs'] = ""
start_DF['Quadrant_Enriched_GOs'] = start_DF['Quadrant_Enriched_GOs'].fillna('')
for idx, row in start_DF.iterrows():
    if pd.isna(row['quadrant']) or pd.isna(row['Athaliana_cleaned_ID']): continue
    quad = row['quadrant']; ath_gene = row['Athaliana_cleaned_ID']
    go_terms_str = row['Quadrant_Enriched_GOs']
    if go_terms_str and isinstance(go_terms_str, str):
        for term_str in go_terms_str.split('; '):
            # Use the new parsing function
            go_id, go_name, p_value = parse_go_term_string_with_pval(term_str)
            if go_id:
                # Store GO info with p-value if available
                go_info_per_quadrant[quad].append((go_id, go_name, p_value))
                go_term_to_genes_per_quad[quad][go_id].add(ath_gene)
print("Finished processing GO terms and p-values.")

# --- 6. Build and Plot Restricted Networks Per Quadrant ---
print("\n--- Building and Plotting Restricted Networks ---")
print(f"Node filter: Including GO terms with >= {MIN_GENES_FOR_NODE} genes AND p_corr <= {GO_PVAL_THRESHOLD}.")
colors = plt.cm.get_cmap('turbo', 25); cluster_colors = {i: colors(i / 25.0) for i in range(1, 26)}; cluster_colors[0] = '#B0B0B0'

if not go_info_per_quadrant: print("No GO terms found. Cannot build networks.")
else:
    for quadrant in sorted(go_info_per_quadrant.keys()):
        go_term_list = go_info_per_quadrant[quadrant]
        if not go_term_list: continue

        print(f"\nProcessing Restricted Network for Quadrant: {quadrant}")
        G_anem = nx.Graph(); G_asag = nx.Graph(); node_data = {}
        valid_nodes_in_quadrant = []
        go_term_data_for_table = [] # For node table output

        # Add Nodes (Applying Stricter Restrictions)
        print(f"  Adding GO term nodes (Min Genes={MIN_GENES_FOR_NODE}, Pval<={GO_PVAL_THRESHOLD})...")
        # Make GO terms unique before processing
        unique_go_terms = {item[0]: item for item in go_term_list}.values()

        for go_id, go_name, p_value in unique_go_terms:
            genes_for_go = go_term_to_genes_per_quad.get(quadrant, {}).get(go_id, set())
            node_gene_count = len(genes_for_go)

            # Apply P-value filter (allow node if p-value missing but gene count sufficient)
            pval_ok = (p_value is not None and p_value <= GO_PVAL_THRESHOLD) or (p_value is None)
            # Apply Gene count filter
            genecount_ok = node_gene_count >= MIN_GENES_FOR_NODE

            if pval_ok and genecount_ok:
                 # Use -log10(p_value) for potential size scaling later, handle p=0 or None
                 if p_value is not None and p_value > 0:
                      significance_score = -math.log10(p_value)
                 else:
                      significance_score = 0 # Assign 0 if p=0 or missing

                 # Scale size by gene count
                 scaled_node_size = node_gene_count * 20 + 50
                 node_label = f"{go_id}\n({go_name})" if go_name != go_id else go_id

                 node_data[go_id] = {'name': node_label, 'size': scaled_node_size,
                                    'genes': genes_for_go, 'p_value': p_value,
                                    'significance': significance_score}
                 G_anem.add_node(go_id, name=node_label, size=scaled_node_size, significance=significance_score)
                 G_asag.add_node(go_id, name=node_label, size=scaled_node_size, significance=significance_score)
                 valid_nodes_in_quadrant.append(go_id)
                 go_term_data_for_table.append({
                     'GO_ID': go_id, 'Name': go_name if go_name != go_id else '',
                     'P_Value': p_value, 'Significance': significance_score,
                     'NodeSize_Scaled': scaled_node_size, 'GeneCount': node_gene_count,
                     'Associated_Ath_Genes': ';'.join(sorted(list(genes_for_go)))
                 })

        print(f"    Added {len(valid_nodes_in_quadrant)} nodes passing filters.")
        if not valid_nodes_in_quadrant: print(f"    No nodes pass filters. Skipping edges/plot."); continue

        # Add Edges (Shared Non-Zero Cluster)
        print("  Adding edges based on shared NON-ZERO motif clusters...")
        edge_data_anem = []; edge_data_asag = []
        can_add_anem_edges = len(anem_gene_to_motifclust) > 0
        can_add_asag_edges = len(asag_gene_to_motifclust) > 0
        node_pairs = []
        for i in range(len(valid_nodes_in_quadrant)):
            for j in range(i + 1, len(valid_nodes_in_quadrant)): node_pairs.append((valid_nodes_in_quadrant[i], valid_nodes_in_quadrant[j]))

        for go_id_a, go_id_b in node_pairs:
            genes_a = node_data[go_id_a]['genes']; genes_b = node_data[go_id_b]['genes']
            # Anem Edges
            if can_add_anem_edges:
                ortho_a = set().union(*(ath_to_anem.get(g, set()) for g in genes_a)); ortho_b = set().union(*(ath_to_anem.get(g, set()) for g in genes_b))
                clusters_a = set(mc[1] for g in ortho_a for mc in anem_gene_to_motifclust.get(g, set())); clusters_b = set(mc[1] for g in ortho_b for mc in anem_gene_to_motifclust.get(g, set()))
                shared_nonzero = {c for c in clusters_a.intersection(clusters_b) if c > 0}
                if shared_nonzero:
                    rep_c = min(shared_nonzero); color = cluster_colors.get(rep_c, 'black')
                    G_anem.add_edge(go_id_a, go_id_b, cluster=rep_c, color=color); edge_data_anem.append({'GO_A': go_id_a, 'GO_B': go_id_b, 'RepresentativeCluster': rep_c, 'SharedClusters': sorted(list(shared_nonzero))})
            # Asag Edges
            if can_add_asag_edges:
                ortho_a = set().union(*(ath_to_asag.get(g, set()) for g in genes_a)); ortho_b = set().union(*(ath_to_asag.get(g, set()) for g in genes_b))
                clusters_a = set(mc[1] for g in ortho_a for mc in asag_gene_to_motifclust.get(g, set())); clusters_b = set(mc[1] for g in ortho_b for mc in asag_gene_to_motifclust.get(g, set()))
                shared_nonzero = {c for c in clusters_a.intersection(clusters_b) if c > 0}
                if shared_nonzero:
                    rep_c = min(shared_nonzero); color = cluster_colors.get(rep_c, 'black')
                    G_asag.add_edge(go_id_a, go_id_b, cluster=rep_c, color=color); edge_data_asag.append({'GO_A': go_id_a, 'GO_B': go_id_b, 'RepresentativeCluster': rep_c, 'SharedClusters': sorted(list(shared_nonzero))})

        print(f"    Added {G_anem.number_of_edges()} Anem edges.")
        print(f"    Added {G_asag.number_of_edges()} Asag edges.")

        # Plotting (Kamada-Kawai with Fallback)
        print("  Plotting restricted networks (Kamada-Kawai)...")
        if not G_anem.nodes() and not G_asag.nodes(): print("    Skipping plot: No nodes."); continue
        fig, axes = plt.subplots(1, 2, figsize=(28, 14))
        fig.suptitle(f'Quadrant {quadrant}: GO Network (>{MIN_GENES_FOR_NODE-1} genes/GO, p<{GO_PVAL_THRESHOLD}) linked by Non-Zero Motif Clusters', fontsize=18)

        # Plot Anem (Left)
        ax1 = axes[0]; ax1.set_title(f'A. nemorensis Links ({G_anem.number_of_edges()} edges)', fontsize=14)
        if G_anem.nodes():
            try: pos_anem = nx.kamada_kawai_layout(G_anem) if nx.is_connected(G_anem) else nx.spring_layout(G_anem, seed=42, k=1.2, iterations=80) # Increase K for disconnected
            except Exception: pos_anem = nx.spring_layout(G_anem, seed=42, k=1.2, iterations=80)
            sizes = [G_anem.nodes[n]['size'] for n in G_anem.nodes()]; labels = {n: G_anem.nodes[n]['name'] for n in G_anem.nodes()}
            colors_e = [G_anem.edges[e]['color'] for e in G_anem.edges()] if G_anem.edges() else []
            nx.draw_networkx_nodes(G_anem, pos_anem, node_size=sizes, node_color='skyblue', alpha=0.8, ax=ax1)
            if G_anem.edges(): nx.draw_networkx_edges(G_anem, pos_anem, edge_color=colors_e, alpha=0.5, width=1.5, ax=ax1)
            nx.draw_networkx_labels(G_anem, pos_anem, labels=labels, font_size=7, ax=ax1) # Slightly smaller font
        ax1.axis('off')

        # Plot Asag (Right)
        ax2 = axes[1]; ax2.set_title(f'A. sagittata Links ({G_asag.number_of_edges()} edges)', fontsize=14)
        if G_asag.nodes():
            try: pos_asag = nx.kamada_kawai_layout(G_asag) if nx.is_connected(G_asag) else nx.spring_layout(G_asag, seed=42, k=1.2, iterations=80)
            except Exception: pos_asag = nx.spring_layout(G_asag, seed=42, k=1.2, iterations=80)
            sizes = [G_asag.nodes[n]['size'] for n in G_asag.nodes()]; labels = {n: G_asag.nodes[n]['name'] for n in G_asag.nodes()}
            colors_e = [G_asag.edges[e]['color'] for e in G_asag.edges()] if G_asag.edges() else []
            nx.draw_networkx_nodes(G_asag, pos_asag, node_size=sizes, node_color='lightcoral', alpha=0.8, ax=ax2)
            if G_asag.edges(): nx.draw_networkx_edges(G_asag, pos_asag, edge_color=colors_e, alpha=0.5, width=1.5, ax=ax2)
            nx.draw_networkx_labels(G_asag, pos_asag, labels=labels, font_size=7, ax=ax1) # Use ax1 font size
        ax2.axis('off')

        # Legend
        present_clusters = set(d['cluster'] for u,v,d in G_anem.edges(data=True)) | set(d['cluster'] for u,v,d in G_asag.edges(data=True)); present_clusters.add(0)
        handles = [plt.Rectangle((0,0),1,1, color=cluster_colors[i]) for i in sorted(list(present_clusters))]
        labels = [f'Cluster {i}' if i > 0 else 'Unassigned' for i in sorted(list(present_clusters))]
        if handles: fig.legend(handles, labels, title="Shared Motif Cluster", loc='lower center', ncol=min(10, len(handles)), fontsize=9) # Adjusted legend

        # Save
        plot_filename = os.path.join(output_dir, f'Quadrant_{quadrant}_Network_KK_Strict.png') # Add Strict
        plt.subplots_adjust(bottom=0.08); plt.savefig(plot_filename, dpi=200, bbox_inches='tight'); plt.close(fig)
        print(f"    Saved restricted network plot: {plot_filename}")

        # Output Restricted Node/Edge Tables
        print("  Saving restricted node and edge tables...")
        if go_term_data_for_table: pd.DataFrame(go_term_data_for_table).to_csv(os.path.join(output_dir, f'Quadrant_{quadrant}_Nodes_Strict.tsv'), sep='\t', index=False)
        if edge_data_anem: edge_df=pd.DataFrame(edge_data_anem); edge_df['SharedClusters']=edge_df['SharedClusters'].apply(lambda x: ';'.join(map(str, x))); edge_df.to_csv(os.path.join(output_dir, f'Quadrant_{quadrant}_Edges_Anem_Strict.tsv'), sep='\t', index=False)
        if edge_data_asag: edge_df=pd.DataFrame(edge_data_asag); edge_df['SharedClusters']=edge_df['SharedClusters'].apply(lambda x: ';'.join(map(str, x))); edge_df.to_csv(os.path.join(output_dir, f'Quadrant_{quadrant}_Edges_Asag_Strict.tsv'), sep='\t', index=False)

print("\n--- Analysis finished ---")
print(f"Script finished at: {pd.Timestamp.now()}")
