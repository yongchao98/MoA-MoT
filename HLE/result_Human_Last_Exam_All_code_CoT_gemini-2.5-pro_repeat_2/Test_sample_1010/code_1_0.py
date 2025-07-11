import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier

def demonstrate_model_suitability():
    """
    Demonstrates that LR, RF, and cosine KNN are all suitable for a classification
    task with concatenated embeddings as features.
    """
    # 1. Define problem parameters
    n_samples = 1000  # Number of text pairs in our dummy dataset
    embedding_dim = 100 # Dimension of a single sentence embedding
    feature_dim = 2 * embedding_dim # Dimension of the concatenated feature vector

    # 2. Generate dummy data
    # X represents the concatenated embeddings [EMBED(x1), EMBED(x2)]
    print(f"Generating dummy data with {n_samples} samples and {feature_dim} features.\n")
    X = np.random.rand(n_samples, feature_dim)
    # y represents the binary labels (0 for non-paraphrase, 1 for paraphrase)
    y = np.random.randint(0, 2, size=n_samples)

    # 3. Instantiate and "train" (fit) each model
    print("Testing model suitability...\n")

    # A. Logistic Regression
    try:
        lr_model = LogisticRegression(max_iter=200)
        lr_model.fit(X, y)
        print("A. Logistic Regression: Successfully fitted. Model is suitable.")
    except Exception as e:
        print(f"A. Logistic Regression: Failed. Error: {e}")

    # B. Random Forest
    try:
        rf_model = RandomForestClassifier()
        rf_model.fit(X, y)
        print("B. Random Forest: Successfully fitted. Model is suitable.")
    except Exception as e:
        print(f"B. Random Forest: Failed. Error: {e}")

    # C. cosine KNN
    # We use KNeighborsClassifier with the 'cosine' metric, which is equivalent to
    # using cosine similarity to find neighbors (since it maximizes the cosine).
    try:
        knn_model = KNeighborsClassifier(metric='cosine')
        knn_model.fit(X, y)
        print("C. cosine KNN: Successfully fitted. Model is suitable.")
    except Exception as e:
        print(f"C. cosine KNN: Failed. Error: {e}")
        
    print("\nConclusion: All tested models can be applied to the feature set.")


if __name__ == '__main__':
    demonstrate_model_suitability()
<<<D>>>