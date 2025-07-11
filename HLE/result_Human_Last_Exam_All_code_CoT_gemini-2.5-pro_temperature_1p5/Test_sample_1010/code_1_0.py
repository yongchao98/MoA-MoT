import numpy as np
import time
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier

def demonstrate_prediction_scalability():
    """
    Demonstrates the prediction time scalability issue of kNN compared to
    Logistic Regression and Random Forest on large datasets.
    """
    embedding_dim = 100
    feature_dim = 2 * embedding_dim  # Concatenated embeddings

    # We will test on datasets of increasing size
    dataset_sizes = [1000, 10000, 50000]

    print("--- Comparing Model Prediction Times vs. Dataset Size ---\n")
    print(f"{'Dataset Size':<15}{'LR Time (s)':<15}{'RF Time (s)':<15}{'kNN Time (s)':<15}")
    print("-" * 60)

    for n_samples in dataset_sizes:
        # 1. Generate a synthetic dataset of size n_samples
        # The data content doesn't matter, only its size and shape
        print(f"Testing with {n_samples} samples...")
        X_train = np.random.rand(n_samples, feature_dim).astype(np.float32)
        y_train = np.random.randint(0, 2, n_samples)
        
        # A single test point for which we'll predict the class
        X_test = np.random.rand(1, feature_dim).astype(np.float32)

        # 2. Train and time Logistic Regression
        lr_model = LogisticRegression(solver='saga', random_state=42)
        lr_model.fit(X_train, y_train)
        start_time = time.time()
        lr_model.predict(X_test)
        lr_time = time.time() - start_time

        # 3. Train and time Random Forest
        rf_model = RandomForestClassifier(n_estimators=10, n_jobs=-1, random_state=42)
        rf_model.fit(X_train, y_train)
        start_time = time.time()
        rf_model.predict(X_test)
        rf_time = time.time() - start_time

        # 4. "Train" and time kNN
        # For kNN, 'fit' just stores the data. The real work is in 'predict'.
        knn_model = KNeighborsClassifier(n_neighbors=5, metric='cosine', n_jobs=-1)
        knn_model.fit(X_train, y_train)
        start_time = time.time()
        knn_model.predict(X_test)
        knn_time = time.time() - start_time
        
        print(f"{n_samples:<15}{lr_time:<15.6f}{rf_time:<15.6f}{knn_time:<15.6f}\n")

    print("\n--- Conclusion ---")
    print("Logistic Regression (LR) and Random Forest (RF) are 'eager learners'.")
    print("They build a compact model during training. Prediction time is fast and constant, regardless of the training data size.")
    print("\nkNN is a 'lazy learner'. It stores the entire training set.")
    print("Prediction requires comparing the test point to *every* training point.")
    print("As shown, its prediction time grows with the dataset size, making it infeasible for 'massive' datasets.")
    print("\nTherefore, cosine KNN is the model we cannot use for this task.")

demonstrate_prediction_scalability()
<<<C>>>