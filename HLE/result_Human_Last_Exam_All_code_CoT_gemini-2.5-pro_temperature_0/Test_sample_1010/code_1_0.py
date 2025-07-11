import numpy as np

def cosine_similarity(v1, v2):
    """Computes the cosine similarity between two vectors."""
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    if norm_v1 == 0 or norm_v2 == 0:
        return 0
    return dot_product / (norm_v1 * norm_v2)

# --- Simulation Setup ---
# Let d=100 for the sentence embeddings
D_EMBED = 100

# 1. Create a large, non-zero mean vector 'c' to simulate anisotropy.
# This represents a common bias found in real-world embedding spaces.
mean_vector = np.ones(D_EMBED) * 5

# 2. Create specific, centered vectors for four different sentences.
# These represent the unique semantic content of each sentence.
np.random.seed(0) # for reproducibility
sent1_specific = np.random.randn(D_EMBED)
sent2_specific = np.random.randn(D_EMBED)
sent3_specific = np.random.randn(D_EMBED)
sent4_specific = np.random.randn(D_EMBED)

# 3. Create the final "anisotropic" embeddings by adding the mean vector.
emb1 = mean_vector + sent1_specific
emb2 = mean_vector + sent2_specific
emb3 = mean_vector + sent3_specific
emb4 = mean_vector + sent4_specific

# 4. Create two feature vectors for the classifier.
# Pair A: (sent1, sent2)
# Pair B: (sent3, sent4)
# These pairs are constructed from completely different underlying random vectors,
# so they should be considered very different.
feature_vector_A = np.concatenate([emb1, emb2])
feature_vector_B = np.concatenate([emb3, emb4])

# 5. Calculate the cosine similarity between the two concatenated feature vectors.
# In a cosine KNN model, this similarity would be used to determine if Pair A
# and Pair B are "neighbors".
sim_anisotropic = cosine_similarity(feature_vector_A, feature_vector_B)

print("--- Demonstration of Cosine KNN Failure with Anisotropic Embeddings ---")
print(f"Dimension of concatenated feature vectors: {len(feature_vector_A)}")
print(f"Cosine similarity between two feature vectors from different pairs: {sim_anisotropic:.6f}")
print("\nObservation: The similarity is very close to 1.0, even though the pairs were generated from completely different random vectors.")
print("This is because the large mean vector dominates the calculation, making all concatenated vectors point in a similar direction.")
print("A KNN model using this metric would not be able to distinguish between different data points, and thus 'cannot at all learn'.")

# --- For comparison, let's see what happens with centered (isotropic) embeddings ---
feature_vector_A_centered = np.concatenate([sent1_specific, sent2_specific])
feature_vector_B_centered = np.concatenate([sent3_specific, sent4_specific])
sim_centered = cosine_similarity(feature_vector_A_centered, feature_vector_B_centered)
print("\n--- Comparison with Centered (Isotropic) Embeddings ---")
print(f"Cosine similarity with centered embeddings: {sim_centered:.6f}")
print("\nObservation: With centered embeddings, the cosine similarity is close to 0, correctly reflecting that the underlying vectors are unrelated.")
print("In this ideal case, cosine KNN would be a suitable model.")