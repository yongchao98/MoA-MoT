import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

def analyze_mimicry_clusters():
    """
    This function simulates similarity data from the described experiment
    and performs hierarchical clustering to group species.

    Even if the data collection method is flawed, this demonstrates the
    analytical step that a researcher would take.
    """
    # Let's assume we have 5 species of Bombus.
    species = [
        'B. ternarius',    # Often red-banded
        'B. borealis',     # Mostly yellow
        'B. terricola',    # Similar pattern to ternarius, but with yellow/white
        'B. pensylvanicus',# Mostly black
        'B. fervidus'      # Mostly yellow, similar to borealis
    ]
    
    print("Step 1: Simulating Similarity Matrix from 'Undergraduate Rankings'.")
    print("Values range from 0 (not similar) to 1 (very similar).\n")
    # This is a hypothetical similarity matrix.
    # High values mean species are considered highly similar.
    # Note the symmetric nature (similarity of A to B is same as B to A).
    similarity_matrix = np.array([
        [1.0, 0.2, 0.8, 0.1, 0.2],  # B. ternarius
        [0.2, 1.0, 0.3, 0.4, 0.9],  # B. borealis
        [0.8, 0.3, 1.0, 0.1, 0.3],  # B. terricola
        [0.1, 0.4, 0.1, 1.0, 0.5],  # B. pensylvanicus
        [0.2, 0.9, 0.3, 0.5, 1.0]   # B. fervidus
    ])

    # Clustering algorithms work with distances, not similarities.
    # A common conversion is: distance = 1 - similarity.
    print("Step 2: Converting Similarity to Distance.")
    print("The 'equation' for this conversion is: distance = 1.0 - similarity\n")
    distance_matrix = 1 - similarity_matrix
    
    # The linkage function requires a condensed distance matrix (a 1D array).
    condensed_distance = squareform(distance_matrix)
    
    print("Step 3: Performing Hierarchical Clustering.")
    print("This process merges the closest species/clusters step-by-step.\n")
    # Using 'average' linkage: distance between clusters is the average of
    # distances between all pairs of items between the two clusters.
    Z = linkage(condensed_distance, 'average')
    
    print("--- Clustering Results (Linkage Matrix) ---")
    print("Each row represents a merge: [cluster_id_1, cluster_id_2, distance, num_items_in_new_cluster]")
    
    for i, row in enumerate(Z):
        cluster_1 = int(row[0])
        cluster_2 = int(row[1])
        distance = row[2]
        num_items = int(row[3])
        
        # Determine names of merged items/clusters
        name_1 = species[cluster_1] if cluster_1 < len(species) else f"cluster_{cluster_1}"
        name_2 = species[cluster_2] if cluster_2 < len(species) else f"cluster_{cluster_2}"
        new_cluster_id = len(species) + i
        
        print(f"\nMerge Step {i+1}:")
        print(f"  - Merging: '{name_1}' (ID: {cluster_1}) and '{name_2}' (ID: {cluster_2})")
        print(f"  - Merge Distance: {distance:.4f}")
        print(f"  - The new cluster (ID: {new_cluster_id}) now contains {num_items} species.")

analyze_mimicry_clusters()