import numpy as np

def find_most_similar_species():
    """
    Analyzes a hypothetical similarity matrix from student ratings to find the
    two species of bumblebees considered most visually similar.
    """
    # Step 1: Define the species. These are common North American bumblebees
    # known to be involved in mimicry complexes.
    species_names = [
        "B. vosnesenskii", # Black-tailed bumblebee
        "B. bifarius",     # Two-form bumblebee
        "B. fervidus",     # Yellow bumblebee
        "B. pensylvanicus" # American bumblebee
    ]

    # Step 2: Define the hypothetical similarity matrix.
    # Scores are out of 10, representing the average from student ratings.
    # The matrix is symmetric, with high values for visually similar species.
    # For example, B. vosnesenskii and B. pensylvanicus are both often
    # yellow and black, so their similarity score is high. B. bifarius often has
    # extensive red markings, making it less similar to the others in this example.
    similarity_matrix = np.array([
        [10.0, 2.1, 8.5, 9.1],  # B. vosnesenskii
        [2.1, 10.0, 3.4, 2.5],  # B. bifarius
        [8.5, 3.4, 10.0, 7.9],  # B. fervidus
        [9.1, 2.5, 7.9, 10.0]   # B. pensylvanicus
    ])

    # Step 3: Find the maximum similarity score.
    # We set the diagonal to a low value to ensure we only find pairs of different species.
    np.fill_diagonal(similarity_matrix, -1.0)
    
    # Find the row and column index of the maximum value in the matrix.
    max_index = np.unravel_index(np.argmax(similarity_matrix, axis=None), similarity_matrix.shape)
    
    # Step 4: Extract the results.
    species1_index, species2_index = max_index
    max_similarity_score = similarity_matrix[species1_index, species2_index]
    
    species1_name = species_names[species1_index]
    species2_name = species_names[species2_index]

    # Step 5: Print the final result in a clear "equation-like" format.
    # This fulfills the requirement to "output each number in the final equation".
    print("Analysis of Similarity Matrix:")
    print(f"The most visually similar pair of species according to the ratings is:")
    # Final Equation Output:
    print(f"Similarity({species1_name}, {species2_name}) = {max_similarity_score}")

find_most_similar_species()