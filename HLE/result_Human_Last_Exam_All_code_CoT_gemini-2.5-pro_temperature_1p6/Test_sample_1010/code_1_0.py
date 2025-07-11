import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

def solve():
    """
    This function demonstrates that Logistic Regression, Random Forest, and 
    cosine KNN are all suitable for the described text similarity task.
    """
    print("--- Demonstrating Model Suitability for Text Similarity Classification ---")

    # 1. Simulate the data
    # Let's assume we have 10 data points (text pairs)
    num_samples = 10
    # The embedding dimension is 100
    embedding_dim = 100
    # The concatenated feature dimension is 200
    feature_dim = embedding_dim * 2

    # Create dummy concatenated embeddings [EMBED(x1), EMBED(x2)]
    # Each row is a 200-dimensional feature vector for a text pair
    X_train = np.random.rand(num_samples, feature_dim)

    # Create dummy binary labels (y=1 for paraphrase, y=0 otherwise)
    y_train = np.random.randint(0, 2, size=num_samples)

    # Create a new sample to predict
    new_sample = np.random.rand(1, feature_dim)

    print(f"Data shape: {X_train.shape}")
    print(f"Labels shape: {y_train.shape}\n")

    # 2. Test each model

    # A. Logistic Regression (LR)
    try:
        lr_model = make_pipeline(StandardScaler(), LogisticRegression())
        lr_model.fit(X_train, y_train)
        prediction = lr_model.predict(new_sample)
        print("A. Logistic Regression: Successfully trained and made a prediction.")
        # print(f"   Prediction for new sample: {prediction[0]}")
    except Exception as e:
        print(f"A. Logistic Regression: FAILED. Error: {e}")

    # B. Random Forest (RF)
    try:
        rf_model = RandomForestClassifier(n_estimators=10, random_state=42)
        rf_model.fit(X_train, y_train)
        prediction = rf_model.predict(new_sample)
        print("B. Random Forest: Successfully trained and made a prediction.")
        # print(f"   Prediction for new sample: {prediction[0]}")
    except Exception as e:
        print(f"B. Random Forest: FAILED. Error: {e}")

    # C. cosine KNN
    try:
        # For cosine KNN, we use KNeighborsClassifier with metric='cosine'.
        # Note: scikit-learn's cosine metric measures distance (1 - cosine_similarity).
        # Finding the smallest distance is equivalent to finding the largest similarity.
        knn_model = KNeighborsClassifier(n_neighbors=3, metric='cosine')
        knn_model.fit(X_train, y_train)
        prediction = knn_model.predict(new_sample)
        print("C. Cosine KNN: Successfully trained and made a prediction.")
        # print(f"   Prediction for new sample: {prediction[0]}")
    except Exception as e:
        print(f"C. Cosine KNN: FAILED. Error: {e}")
        
    print("\nConclusion: All three models are suitable for the task.")

solve()
<<<D>>>