import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier

def demonstrate_model_suitability():
    """
    This function demonstrates that Logistic Regression, Random Forest, and 
    Cosine k-NN are all suitable for predicting similarity from concatenated embeddings.
    """
    # 1. Define problem parameters
    n_samples = 1000  # Number of text pairs in our dataset
    d = 100           # Dimension of each sentence embedding
    feature_dim = 2 * d # Dimension of the concatenated feature vector

    print(f"Setting up a simulated classification problem:")
    print(f"Number of samples: {n_samples}")
    print(f"Embedding dimension (d): {d}")
    print(f"Final feature vector dimension (2*d): {feature_dim}\n")

    # 2. Generate synthetic data to mimic the problem
    # EMBED(x1): Random embeddings for the first text in each pair
    embeddings1 = np.random.rand(n_samples, d)
    # EMBED(x2): Random embeddings for the second text in each pair
    embeddings2 = np.random.rand(n_samples, d)
    
    # Concatenated features: [EMBED(x1), EMBED(x2)]
    X_features = np.concatenate([embeddings1, embeddings2], axis=1)
    
    # y: Random binary labels (0 for not paraphrase, 1 for paraphrase)
    y_labels = np.random.randint(0, 2, size=n_samples)

    print(f"Generated synthetic data:")
    print(f"Feature matrix shape (X): {X_features.shape}")
    print(f"Label vector shape (y): {y_labels.shape}\n")

    # 3. Instantiate and train each model
    # A single data point to test prediction
    test_sample = X_features[0].reshape(1, -1)

    # Model A: Logistic Regression
    try:
        lr_model = LogisticRegression(max_iter=100)
        lr_model.fit(X_features, y_labels)
        prediction = lr_model.predict(test_sample)
        print(f"A. Logistic Regression: SUITABLE")
        print(f"   - Successfully trained the model.")
        print(f"   - Prediction for a sample point: {prediction[0]}")
    except Exception as e:
        print(f"A. Logistic Regression: FAILED with error {e}")

    # Model B: Random Forest
    try:
        rf_model = RandomForestClassifier(n_estimators=10)
        rf_model.fit(X_features, y_labels)
        prediction = rf_model.predict(test_sample)
        print(f"B. Random Forest: SUITABLE")
        print(f"   - Successfully trained the model.")
        print(f"   - Prediction for a sample point: {prediction[0]}")
    except Exception as e:
        print(f"B. Random Forest: FAILED with error {e}")

    # Model C: cosine k-NN
    try:
        # For k-NN, 'cosine' is a valid metric to measure distance between points.
        knn_model = KNeighborsClassifier(n_neighbors=5, metric='cosine', algorithm='brute')
        knn_model.fit(X_features, y_labels)
        prediction = knn_model.predict(test_sample)
        print(f"C. Cosine k-NN: SUITABLE")
        print(f"   - Successfully trained the model.")
        print(f"   - Prediction for a sample point: {prediction[0]}")
    except Exception as e:
        print(f"C. Cosine k-NN: FAILED with error {e}")
        
    print("\nConclusion: All three models can be instantiated, trained, and used for prediction on the given data structure.")

if __name__ == '__main__':
    demonstrate_model_suitability()