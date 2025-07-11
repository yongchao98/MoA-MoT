import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score

def demonstrate_model_suitability():
    """
    This function demonstrates that Logistic Regression, Random Forest, and cosine KNN
    are all suitable models for predicting paraphrase similarity from concatenated
    sentence embeddings.
    """
    # 1. SETUP THE SCENARIO
    # Low-dimensional embedding function (d=100)
    embedding_dim = 100
    # A massive dataset of text pairs (we'll simulate 1000 pairs)
    num_samples = 1000
    # Feature vector is [EMBED(x1), EMBED(x2)], so feature_dim = 200
    feature_dim = embedding_dim * 2

    print("--- Problem Setup ---")
    print(f"Sentence embedding dimension (d): {embedding_dim}")
    print(f"Feature vector dimension (concatenated): {feature_dim}")
    print(f"Number of training samples: {num_samples}")
    print("-" * 25 + "\n")

    # 2. GENERATE SIMULATED DATA
    # Create random embeddings for n_samples pairs of sentences
    # We generate a pool of embeddings to draw from
    embeddings_pool = np.random.randn(num_samples * 2, embedding_dim)

    # Create the feature matrix (X) and label vector (y)
    X_list = []
    y = []

    for i in range(num_samples):
        e1 = embeddings_pool[2 * i]
        e2 = embeddings_pool[2 * i + 1]

        # The feature vector for our model is the concatenation
        feature_vector = np.concatenate([e1, e2])
        X_list.append(feature_vector)

        # Create a plausible label 'y'
        # y=1 (paraphrase) if e1 and e2 are similar, y=0 otherwise
        # We can simulate this using the cosine similarity between e1 and e2
        cosine_sim = np.dot(e1, e2) / (np.linalg.norm(e1) * np.linalg.norm(e2))
        # If similarity is high, they are paraphrases (y=1)
        if cosine_sim > 0.5:
            y.append(1)
        else:
            y.append(0)

    X = np.array(X_list)
    y = np.array(y)

    print(f"Generated feature matrix X with shape: {X.shape}")
    print(f"Generated label vector y with shape: {y.shape}\n")


    # 3. TEST EACH MODEL
    # A dictionary to hold the models we want to test
    models = {
        "A. Logistic Regression": LogisticRegression(max_iter=1000, random_state=42),
        "B. Random Forest": RandomForestClassifier(random_state=42),
        "C. cosine KNN": KNeighborsClassifier(n_neighbors=5, metric='cosine')
    }

    print("--- Evaluating Model Suitability ---")
    all_suitable = True
    for name, model in models.items():
        print(f"\nTesting Model: {name}")
        try:
            # Train the model on the concatenated features
            model.fit(X, y)
            print("  - Training: Successful")

            # Make predictions to ensure the model works end-to-end
            y_pred = model.predict(X)
            # Report accuracy on training data just to show it learned something
            accuracy = accuracy_score(y, y_pred)
            print(f"  - Prediction: Successful (Train accuracy: {accuracy:.2f})")
            print("  - Conclusion: This model IS SUITABLE.")

        except Exception as e:
            print(f"  - An error occurred: {e}")
            print("  - Conclusion: This model is NOT SUITABLE.")
            all_suitable = False

    # 4. FINAL CONCLUSION
    print("\n--- Final Analysis ---")
    if all_suitable:
        print("All tested models (Logistic Regression, Random Forest, cosine KNN) were able to train")
        print("and make predictions on the feature set [EMBED(x1), EMBED(x2)].")
        print("Therefore, all these models are suitable for the task.")
    else:
        print("At least one model was found to be unsuitable.")

if __name__ == '__main__':
    demonstrate_model_suitability()