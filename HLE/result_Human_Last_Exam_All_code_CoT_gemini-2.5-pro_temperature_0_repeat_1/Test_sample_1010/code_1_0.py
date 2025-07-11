import numpy as np

def cosine_similarity(v1, v2):
    """Computes the cosine similarity between two vectors."""
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    if norm_v1 == 0 or norm_v2 == 0:
        return 0
    return dot_product / (norm_v1 * norm_v2)

def run_demonstration():
    """
    Demonstrates why cosine KNN is unsuitable for this task.
    We define a toy paraphrase rule and show that the cosine similarity
    of concatenated vectors does not reflect the similarity of labels.
    """
    print("--- Toy Example Demonstration ---")
    print("Rule: A pair of vectors (e1, e2) is a paraphrase (Label=1) if the distance ||e1 - e2|| < 1.0, otherwise Label=0.\n")

    # --- Data Point A (Label = 1) ---
    e1a = np.array([1.0, 0.0])
    e2a = np.array([1.0, 0.8])
    dist_a = np.linalg.norm(e1a - e2a)
    label_a = 1 if dist_a < 1.0 else 0
    v_a = np.concatenate([e1a, e2a])
    print(f"Point A: e1a = {e1a}, e2a = {e2a}")
    print(f"||e1a - e2a|| = {dist_a:.2f} -> Label = {label_a}")
    print(f"Concatenated vector v_a = {v_a}\n")

    # --- Data Point B (Label = 1) ---
    # This point also has Label=1, so it should be "close" to Point A.
    e1b = np.array([0.0, 1.0])
    e2b = np.array([0.8, 1.0])
    dist_b = np.linalg.norm(e1b - e2b)
    label_b = 1 if dist_b < 1.0 else 0
    v_b = np.concatenate([e1b, e2b])
    print(f"Point B: e1b = {e1b}, e2b = {e2b}")
    print(f"||e1b - e2b|| = {dist_b:.2f} -> Label = {label_b}")
    print(f"Concatenated vector v_b = {v_b}\n")

    # --- Data Point C (Label = 0) ---
    # This point has Label=0, so it should be "far" from Point A.
    e1c = np.array([1.0, 0.2])
    e2c = np.array([0.0, 1.0])
    dist_c = np.linalg.norm(e1c - e2c)
    label_c = 1 if dist_c < 1.0 else 0
    v_c = np.concatenate([e1c, e2c])
    print(f"Point C: e1c = {e1c}, e2c = {e2c}")
    print(f"||e1c - e2c|| = {dist_c:.2f} -> Label = {label_c}")
    print(f"Concatenated vector v_c = {v_c}\n")

    # --- Calculate Similarities ---
    print("--- Cosine Similarities between Concatenated Vectors ---")
    
    # Similarity between two points with the SAME label (A and B)
    sim_ab = cosine_similarity(v_a, v_b)
    print(f"Similarity(v_a, v_b) where Label(A)=1, Label(B)=1")
    # Final equation for sim(v_a, v_b)
    print(f"Equation: dot({v_a}, {v_b}) / (||{v_a}|| * ||{v_b}||)")
    print(f"Result: {np.dot(v_a, v_b):.2f} / ({np.linalg.norm(v_a):.2f} * {np.linalg.norm(v_b):.2f}) = {sim_ab:.4f}\n")

    # Similarity between two points with DIFFERENT labels (A and C)
    sim_ac = cosine_similarity(v_a, v_c)
    print(f"Similarity(v_a, v_c) where Label(A)=1, Label(C)=0")
    # Final equation for sim(v_a, v_c)
    print(f"Equation: dot({v_a}, {v_c}) / (||{v_a}|| * ||{v_c}||)")
    print(f"Result: {np.dot(v_a, v_c):.2f} / ({np.linalg.norm(v_a):.2f} * {np.linalg.norm(v_c):.2f}) = {sim_ac:.4f}\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    print(f"Point A (Label 1) is more similar to Point C (Label 0) with similarity {sim_ac:.4f}")
    print(f"than it is to Point B (Label 1) with similarity {sim_ab:.4f}.")
    print("A KNN classifier would therefore likely misclassify Point A based on its 'neighbors'.")
    print("This is because the cosine similarity of concatenated vectors is a poor metric for paraphrase detection.")

if __name__ == '__main__':
    run_demonstration()
<<<C>>>