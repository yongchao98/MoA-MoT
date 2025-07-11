import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

def demonstrate_model_suitability():
    """
    This function simulates the paraphrase detection task to demonstrate that
    Logistic Regression, Random Forest, and cosine KNN are all suitable models.
    """
    # 1. Define problem parameters
    embedding_dim = 100
    feature_dim = 2 * embedding_dim  # Concatenated embeddings
    num_samples = 1000

    print(f"Simulating a dataset for paraphrase detection.")
    print(f"Embedding dimension (d): {embedding_dim}")
    print(f"Feature vector dimension (concatenated): {feature_dim}")
    print("-" * 50)

    # 2. Simulate the data
    # X represents the [EMBED(x1), EMBED(x2)] vectors
    X = np.random.rand(num_samples, feature_dim)
    # y represents the labels (1 for paraphrase, 0 otherwise)
    # We generate a simple rule for y to make it learnable
    # e.g., paraphrase if the dot product of the first 10 features is high
    pseudo_similarity = np.sum(X[:, :10] * X[:, embedding_dim:embedding_dim+10], axis=1)
    y = (pseudo_similarity > 2.5).astype(int)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

    # 3. Instantiate the models
    models = {
        "A. Logistic Regression": LogisticRegression(max_iter=200, solver='liblinear'),
        "B. Random Forest": RandomForestClassifier(random_state=42),
        "C. cosine KNN": KNeighborsClassifier(metric='cosine')
    }

    # 4. Train and evaluate each model to prove they are usable
    print("Testing each model's ability to train on the data...\n")
    for name, model in models.items():
        try:
            # Train (fit) the model on the training data
            model.fit(X_train, y_train)

            # Predict on the test data
            y_pred = model.predict(X_test)

            # Calculate the accuracy score (a resulting number)
            score = accuracy_score(y_test, y_pred)

            print(f"Model: {name}")
            print(f"Result: SUCCESS. The model trained and made predictions.")
            print(f"A resulting number (Test Accuracy): {score:.4f}")
            print("-" * 25)

        except Exception as e:
            print(f"Model: {name}")
            print(f"Result: FAILED. This model is not suitable.")
            print(f"Error: {e}")
            print("-" * 25)

if __name__ == '__main__':
    demonstrate_model_suitability()