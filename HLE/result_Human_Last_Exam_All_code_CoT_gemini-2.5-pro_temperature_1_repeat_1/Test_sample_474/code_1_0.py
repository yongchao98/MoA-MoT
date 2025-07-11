import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import itertools

def simulate_bee_similarity_study():
    """
    This script simulates the researcher's experiment of clustering Bombus species
    based on similarity rankings from untrained undergraduates.
    """
    # 1. Define our hypothetical Bombus species
    species = [
        "B. terrestris",  # Buff-tailed (white tail)
        "B. lucorum",     # White-tailed
        "B. lapidarius",  # Red-tailed
        "B. ruderarius",  # Red-shanked carder bee (some red)
        "B. pascuorum",   # Common carder bee (all ginger)
        "B. hortorum"      # Garden bumblebee (long face, white tail)
    ]
    num_species = len(species)
    print(f"Species being studied: {species}\n")

    # 2. Define ground-truth mimicry groups (for simulation purposes).
    # This information is ONLY for creating realistic simulation data.
    # The clustering algorithm will not have access to this truth.
    # Group 1: White-tailed mimics
    # Group 2: Red-tailed mimics
    # Group 3: All-ginger mimic
    true_groups = {
        "B. terrestris": 1, "B. lucorum": 1, "B. hortorum": 1,
        "B. lapidarius": 2, "B. ruderarius": 2,
        "B. pascuorum": 3
    }

    # 3. Simulate ratings from 20 undergraduates
    num_undergrads = 20
    # Create a similarity matrix to store the sum of scores.
    # Similarity is rated on a scale of 1 (very different) to 7 (very similar).
    similarity_sum = np.zeros((num_species, num_species))

    # Iterate through all unique pairs of species
    for i, j in itertools.combinations(range(num_species), 2):
        sp1, sp2 = species[i], species[j]
        # Base similarity on the 'true' groups to make the simulation realistic.
        if true_groups[sp1] == true_groups[sp2]:
            # If in the same mimicry group, ratings should be high
            base_similarity = 6
        else:
            # If in different groups, ratings should be low
            base_similarity = 2

        # Add random noise for each of the 20 undergrads' scores
        student_scores = base_similarity + np.random.randint(-1, 2, size=num_undergrads)
        # Ensure scores stay within the 1-7 range
        student_scores = np.clip(student_scores, 1, 7)
        total_score = np.sum(student_scores)

        similarity_sum[i, j] = total_score
        similarity_sum[j, i] = total_score # The matrix is symmetric

    # 4. Calculate the average similarity and convert it to a distance matrix
    # Hierarchical clustering requires a distance matrix, where high values mean dissimilar.
    avg_similarity = similarity_sum / num_undergrads
    # Convert similarity to distance (e.g., distance = max_similarity - similarity)
    distance_matrix = 7 - avg_similarity
    np.fill_diagonal(distance_matrix, 0)

    # Convert the square distance matrix to a condensed distance matrix for scipy
    condensed_distance = squareform(distance_matrix)

    print("Average Similarity Matrix (from 1 to 7, higher is more similar):")
    print(np.round(avg_similarity, 2))
    print("\nConverted Distance Matrix (lower is more similar):")
    print(np.round(distance_matrix, 2))
    print("\n")

    # 5. Perform Hierarchical Clustering on the distance matrix
    # The 'average' method uses the average distance between elements of each cluster.
    linked = linkage(condensed_distance, 'average')

    # 6. Form flat clusters from the hierarchy using a distance threshold.
    # A low threshold means only very similar items will be grouped.
    # We'll set a threshold of 3.5 (meaning avg similarity must be > 3.5 to cluster).
    cluster_threshold = 3.5
    clusters = fcluster(linked, cluster_threshold, criterion='distance')

    # 7. Print the final clusters
    print(f"Final Clusters (based on a distance threshold of {cluster_threshold}):")
    result = {}
    for i, species_name in enumerate(species):
        cluster_id = clusters[i]
        if cluster_id not in result:
            result[cluster_id] = []
        result[cluster_id].append(species_name)

    for cluster_id, members in sorted(result.items()):
        print(f"Cluster {cluster_id}: {members}")

# Run the simulation of the experiment
simulate_bee_similarity_study()