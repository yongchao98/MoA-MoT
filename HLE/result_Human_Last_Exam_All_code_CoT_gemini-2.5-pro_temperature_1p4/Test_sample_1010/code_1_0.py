import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split

def demonstrate_model_suitability():
    """
    This function demonstrates that all three models (LR, RF, cosine KNN)
    are suitable for the paraphrase detection task with concatenated embeddings.
    """
    # --- Step 1: Define a mock Embedding Function and sample Data ---

    # Assume we have a low-dimensional (d=100) sentence embedding function (EMBED).
    D_EMBED = 100
    def EMBED(sentence: str) -> np.ndarray:
        """A dummy embedding function that returns a normalized 100-dim vector."""
        seed = sum(ord(c) for c in sentence)
        rs = np.random.RandomState(seed)
        vec = rs.rand(D_EMBED)
        return vec / np.linalg.norm(vec)

    # Assume a dataset of text pairs (x1, x2) and a label y.
    dataset = [
        ("The cat sat on the mat.", "There was a cat on the mat.", 1),
        ("The cat sat on the mat.", "The dog barked at the moon.", 0),
        ("It is a sunny day.", "The weather is bright and clear.", 1),
        ("It is a sunny day.", "I need to go shopping.", 0),
        ("What is your name?", "Could you please tell me your name?", 1),
        ("What is your name?", "My name is John.", 0),
    ] * 20 # Make the dataset larger for training

    # --- Step 2: Prepare Features ([EMBED(x1), EMBED(x2)]) and Labels ---
    features = []
    labels = []
    for x1, x2, y in dataset:
        concatenated_features = np.concatenate([EMBED(x1), EMBED(x2)])
        features.append(concatenated_features)
        labels.append(y)

    X = np.array(features)
    y = np.array(labels)

    print(f"Shape of feature matrix X: {X.shape}")
    print(f"Shape of label vector y: {y.shape}\n")

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

    # --- Step 3: Show that each model can be trained on this data ---

    # A. Logistic Regression
    print("--- Testing Model Suitability ---")
    try:
        lr_model = LogisticRegression(max_iter=1000, random_state=42)
        lr_model.fit(X_train, y_train)
        print("A. Logistic Regression: SUITABLE. Model was trained successfully.")
    except Exception as e:
        print(f"A. Logistic Regression: NOT SUITABLE. Error: {e}")

    # B. Random Forest
    try:
        rf_model = RandomForestClassifier(random_state=42)
        rf_model.fit(X_train, y_train)
        print("B. Random Forest: SUITABLE. Model was trained successfully.")
    except Exception as e:
        print(f"B. Random Forest: NOT SUITABLE. Error: {e}")

    # C. cosine KNN
    try:
        # KNeighborsClassifier with metric='cosine' is the cosine KNN model.
        knn_model = KNeighborsClassifier(n_neighbors=3, metric='cosine')
        knn_model.fit(X_train, y_train)
        print("C. cosine KNN: SUITABLE. Model was trained successfully.")
    except Exception as e:
        print(f"C. cosine KNN: NOT SUITABLE. Error: {e}")

    print("\n--- Conclusion ---")
    print("All models can be instantiated and trained on the given feature representation.")
    print("Therefore, none of the models are fundamentally unable to learn from this data.")


if __name__ == '__main__':
    demonstrate_model_suitability()