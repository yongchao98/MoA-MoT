import numpy as np
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics.pairwise import cosine_similarity

# Plan:
# 1. Simulate a simple 2D embedding space for words like 'king', 'monarch', 'queen'.
# 2. Create training pairs where two pairs are very "close" in the concatenated feature space
#    but have opposite paraphrase labels.
# 3. Show that because of this, a KNN classifier gets confused and makes unreliable predictions,
#    as it cannot learn the actual rule for what makes a paraphrase.

# 1. Simulate a simple 2D embedding space (d=2 for easy demonstration)
# Imagine these are the outputs of our EMBED() function.
embeddings = {
    "king": np.array([0.9, 0.1]),
    "monarch": np.array([0.85, 0.15]), # Semantically similar to 'king'
    "queen": np.array([0.8, 0.2]),     # Also semantically similar to 'king'
    "car": np.array([-0.5, 0.8]),
    "auto": np.array([-0.45, 0.85])  # Semantically similar to 'car'
}

# Define our paraphrase rule: y=1 if cosine similarity > 0.99, else y=0
def is_paraphrase(emb1, emb2):
    similarity = cosine_similarity(emb1.reshape(1, -1), emb2.reshape(1, -1))[0, 0]
    return 1 if similarity > 0.99 else 0

# 2. Create a training set that demonstrates the problem for KNN
# Training Pair A: A paraphrase about royalty
x1_A = embeddings["king"]
x2_A = embeddings["monarch"]
y_A = is_paraphrase(x1_A, x2_A) # Label is 1
X_A = np.concatenate([x1_A, x2_A]) # Concatenated feature vector

# Training Pair B: A non-paraphrase, but topically similar to Pair A
x1_B = embeddings["king"]
x2_B = embeddings["queen"]
y_B = is_paraphrase(x1_B, x2_B) # Label is 0
X_B = np.concatenate([x1_B, x2_B])

# Training Pair C: A paraphrase about cars, topically different from A and B
x1_C = embeddings["car"]
x2_C = embeddings["auto"]
y_C = is_paraphrase(x1_C, x2_C) # Label is 1
X_C = np.concatenate([x1_C, x2_C])

# Our full training data
X_train = np.array([X_A, X_B, X_C])
y_train = np.array([y_A, y_B, y_C])

# 3. Show that two points with different labels are "close" in the feature space
# We will use cosine similarity, as the model is "cosine KNN". High similarity means low distance.
sim_A_B = cosine_similarity(X_A.reshape(1, -1), X_B.reshape(1, -1))[0, 0]
sim_A_C = cosine_similarity(X_A.reshape(1, -1), X_C.reshape(1, -1))[0, 0]

print("--- The Core Problem for KNN ---")
print(f"Training Pair A: ['king', 'monarch'] -> Concatenated Vector, Label: {y_A} (Paraphrase)")
print(f"Training Pair B: ['king', 'queen']   -> Concatenated Vector, Label: {y_B} (Not Paraphrase)")
print(f"Training Pair C: ['car', 'auto']     -> Concatenated Vector, Label: {y_C} (Paraphrase)")
print("-" * 30)
print("Let's measure the 'nearness' (cosine similarity) between these training vectors:")
print(f"Similarity(Vector A, Vector B): {sim_A_B:.4f}")
print(f"Similarity(Vector A, Vector C): {sim_A_C:.4f}")
print("-" * 30)
print("Observation:")
print("Vector A and Vector B are extremely similar (neighbors) in the feature space.")
print("However, they have OPPOSITE labels (1 and 0).")
print("This violates the fundamental assumption of KNN, which expects neighbors to have the same label.")
print("\n--- KNN Prediction Demonstration ---")
# A test point that is a paraphrase and very similar to Pair A
x1_test = embeddings["king"]
x2_test = np.array([0.86, 0.14]) # A new word, very similar to 'monarch'
y_test_true = is_paraphrase(x1_test, x2_test)
X_test = np.concatenate([x1_test, x2_test])

# Because the test point is so similar to both A and B, the 1-NN prediction is unstable.
# In this specific case, let's see which is closer.
sim_test_A = cosine_similarity(X_test.reshape(1, -1), X_A.reshape(1, -1))[0, 0]
sim_test_B = cosine_similarity(X_test.reshape(1, -1), X_B.reshape(1, -1))[0, 0]

print(f"Test point's true label is: {y_test_true}")
print(f"Similarity of Test point to A (Label 1): {sim_test_A:.4f}")
print(f"Similarity of Test point to B (Label 0): {sim_test_B:.4f}")

# The prediction will be the label of the nearest neighbor
neighbor_label = y_A if sim_test_A > sim_test_B else y_B
print(f"The nearest neighbor is {'A' if sim_test_A > sim_test_B else 'B'}, so the predicted label is {neighbor_label}.")

print("\nFinal Conclusion:")
print("The model's prediction depends on tiny variations in semantic similarity, not on learning the actual rule of paraphrasing.")
print("It is not 'learning' a function to compare the two parts of the vector; it is only looking at the similarity of entire concatenated vectors.")
print("Therefore, cosine KNN is not at all suitable for this task with this feature representation.")
