import numpy as np

def run_demonstration():
    """
    This function demonstrates why 'cosine KNN' on concatenated embeddings is
    unsuitable for paraphrase detection.
    """
    # Use a small dimension for embeddings for clarity. The logic holds for d=100.
    # We create mock concept vectors that are normalized (length = 1).
    def normalize(v):
        return v / np.linalg.norm(v)

    # Concepts: 'cat' is similar to 'feline', 'dog' is similar to 'canine'.
    # 'cat' and 'dog' topics are dissimilar (orthogonal).
    V_cat = normalize(np.array([1.0, 0.1]))
    V_feline = normalize(np.array([1.0, 0.2]))
    V_dog = normalize(np.array([0.1, 1.0]))
    V_canine = normalize(np.array([0.2, 1.0]))

    # Mock EMBED function creates a sentence embedding from concepts.
    def EMBED(text):
        if "cat" in text:
            return V_cat
        if "feline" in text:
            return V_feline
        if "dog" in text:
            return V_dog
        if "canine" in text:
            return V_canine
        return np.zeros_like(V_cat)

    # --- Define three pairs and their feature vectors ---

    # Pair 1: Paraphrase (y=1) on the 'cat' topic.
    E1a = EMBED("the cat")
    E1b = EMBED("the feline")
    F1 = np.concatenate([E1a, E1b])

    # Pair 2: Paraphrase (y=1) on the 'dog' topic.
    E2a = EMBED("the dog")
    E2b = EMBED("the canine")
    F2 = np.concatenate([E2a, E2b])

    # Pair 3: Not a paraphrase (y=0), mixing 'cat' and 'dog' topics.
    # It shares its first sentence with Pair 1.
    E3a = EMBED("the cat")
    E3b = EMBED("the dog")
    F3 = np.concatenate([E3a, E3b])
    
    # --- Define distance calculation ---
    def cosine_similarity(v1, v2):
        return np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))

    def cosine_distance(v1, v2):
        return 1 - cosine_similarity(v1, v2)
        
    # --- The Argument ---
    # For KNN to work, points with the same label should be close.
    # Let's check if F1 (y=1) is closer to F2 (y=1) or F3 (y=0).
    dist_1_to_2 = cosine_distance(F1, F2)
    dist_1_to_3 = cosine_distance(F1, F3)

    print("This script shows why 'cosine KNN' is unsuitable for this paraphrase task.")
    print("The distance metric for KNN on concatenated vectors is misaligned with the task goal.\n")
    print("Let F1 = [E('cat'), E('feline')] --- Paraphrase, y=1")
    print("Let F2 = [E('dog'), E('canine')] --- Paraphrase, y=1")
    print("Let F3 = [E('cat'), E('dog')] --- Not Paraphrase, y=0\n")
    print("For KNN to work, distance(F1, F2) should be smaller than distance(F1, F3).\n")

    print("--- Calculating Distance(F1, F2) ---")
    sim_1_2 = cosine_similarity(F1, F2)
    # The norm of a concatenation of two unit vectors is sqrt(1^2 + 1^2) = sqrt(2)
    norm_F = np.sqrt(2)
    dot_E1a_E2a = np.dot(E1a, E2a)
    dot_E1b_E2b = np.dot(E1b, E2b)
    print(f"Similarity = (dot(E1a, E2a) + dot(E1b, E2b)) / (norm(F1) * norm(F2))")
    print(f"           = ({dot_E1a_E2a:.2f} + {dot_E1b_E2b:.2f}) / ({np.linalg.norm(F1):.2f} * {np.linalg.norm(F2):.2f}) = {sim_1_2:.4f}")
    print(f"Distance(F1, F2) = 1.0 - {sim_1_2:.4f} = {dist_1_to_2:.4f}\n")

    print("--- Calculating Distance(F1, F3) ---")
    sim_1_3 = cosine_similarity(F1, F3)
    dot_E1a_E3a = np.dot(E1a, E3a)
    dot_E1b_E3b = np.dot(E1b, E3b)
    print(f"Similarity = (dot(E1a, E3a) + dot(E1b, E3b)) / (norm(F1) * norm(F3))")
    print(f"           = ({dot_E1a_E3a:.2f} + {dot_E1b_E3b:.2f}) / ({np.linalg.norm(F1):.2f} * {np.linalg.norm(F3):.2f}) = {sim_1_3:.4f}")
    print(f"Distance(F1, F3) = 1.0 - {sim_1_3:.4f} = {dist_1_to_3:.4f}\n")

    print("--- Conclusion ---")
    print(f"Distance between two paraphrase pairs (F1, F2): {dist_1_to_2:.4f}")
    print(f"Distance between a paraphrase and non-paraphrase pair (F1, F3): {dist_1_to_3:.4f}\n")
    
    if dist_1_to_2 > dist_1_to_3:
        print("Result: The two paraphrase pairs are MORE distant than the paraphrase/non-paraphrase pair.")
        print("A KNN model would consider F3 a closer neighbor to F1 than F2 is.")
        print("This is because the distance metric rewards pairs that are topically similar, not pairs that exhibit the internal property of being paraphrases.")
        print("Therefore, the model is not suitable for learning this task.")
    else:
        # This case shouldn't happen with the chosen vectors, but is included for completeness.
        print("While the distances in this run didn't show the issue, the underlying logic holds:")
        print("The distance metric is fundamentally misaligned with the classification objective.")

if __name__ == "__main__":
    run_demonstration()