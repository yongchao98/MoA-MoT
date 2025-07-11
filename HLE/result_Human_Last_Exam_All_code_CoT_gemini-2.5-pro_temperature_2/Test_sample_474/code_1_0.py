import numpy as np

def evaluate_mimicry_method():
    """
    This function demonstrates the flaw in using human perception to assess
    mimicry syndromes that are targeted at avian predators.
    """
    print("Evaluating a research method for Bombus mimicry clustering.")
    print("----------------------------------------------------------\n")
    print("The proposed method's validity depends on whether human vision is a good proxy for the vision of the primary predator, typically birds.")
    print("Birds have different visual systems (often tetrachromatic, seeing UV light), so human perception can be misleading.\n")

    # Define our species
    species = ['Bombus A (Highland)', 'Bombus B (Highland)', 'Bombus C (Lowland)']

    # 1. The "True" Ecological Similarity (as seen by a bird predator)
    # Let's assume A and B are true Müllerian mimics, but C is not.
    true_similarity = np.array([
        [1.0, 0.9, 0.1],  # A is identical to A, very similar to B, not similar to C
        [0.9, 1.0, 0.15], # B is very similar to A, identical to B, not similar to C
        [0.1, 0.15, 1.0]  # C is not similar to A or B
    ])

    print("Step 1: Define a 'True' Similarity Matrix (Bird's Perception)")
    print("This represents the actual ecological similarities that drive mimicry.")
    print(species)
    print(true_similarity)
    print("\n")

    # 2. The "Human-Perceived" Similarity (undergraduate ranking)
    # We simulate this by adding 'noise' to the true matrix. This noise
    # represents the error introduced by using the wrong visual system and untrained observers.
    # We will seed the random number generator for reproducible results.
    np.random.seed(42)
    # The noise will be a random value between -0.2 and 0.2
    noise = (np.random.rand(3, 3) - 0.5) * 0.4
    human_perceived_similarity = true_similarity + noise
    # Clamp values to the valid [0, 1] range for similarity
    human_perceived_similarity = np.clip(human_perceived_similarity, 0, 1)
    # Ensure the diagonal is 1.0 (a species is always identical to itself)
    np.fill_diagonal(human_perceived_similarity, 1.0)


    print("Step 2: Simulate a 'Human-Perceived' Similarity Matrix")
    print("This matrix includes errors because human vision differs from bird vision.")
    print(species)
    print(np.round(human_perceived_similarity, 4))
    print("\n")

    # 3. Quantify the error
    # We use Mean Squared Error (MSE) to measure the difference.
    squared_errors = (true_similarity - human_perceived_similarity) ** 2
    mse = np.mean(squared_errors)

    print("Step 3: Calculate the error between the two matrices.")
    print("A high error rate shows that the human-based method is inaccurate.\n")

    # Show a sample calculation for the final equation as requested.
    # Let's show the error for the A-B similarity.
    true_ab = true_similarity[0, 1]
    human_ab = human_perceived_similarity[0, 1]
    error_ab_sq = squared_errors[0, 1]
    
    print("Example Calculation (Error for Species A-B Similarity):")
    # I am printing each number of the equation as requested
    print(f"True Similarity: {true_ab}")
    print(f"Human-Perceived Similarity: {human_ab:.4f}")
    print(f"Equation: (true_value - perceived_value)^2 = error")
    print(f"Calculation: ({true_ab} - {human_ab:.4f})^2 = {error_ab_sq:.4f}")
    print("\n")
    
    print(f"The overall Mean Squared Error (MSE) between the true and human-perceived matrices is: {mse:.4f}")

    print("\n----------------------------------------------------------")
    print("Conclusion: Because the error is significant, the human-based method is not ecologically valid.")

if __name__ == '__main__':
    evaluate_mimicry_method()