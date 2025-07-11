import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier

def demonstrate_model_suitability():
    """
    This function demonstrates that LR, RF, and Cosine KNN can all be used
    for a binary classification task with concatenated embeddings as features.
    """
    # d: dimension of a single sentence embedding
    d = 100
    
    # n_samples: number of text pairs in our dummy training set
    n_samples = 50
    
    # 1. Create a dummy dataset
    # X_train represents the concatenated embeddings [EMBED(x1), EMBED(x2)]
    # Shape is (n_samples, 2 * d)
    print(f"Creating a dummy dataset with {n_samples} samples.")
    X_train = np.random.rand(n_samples, 2 * d)
    
    # y_train represents the binary labels (0 or 1)
    # Shape is (n_samples,)
    y_train = np.random.randint(0, 2, n_samples)
    
    print(f"Input feature vector shape for training: {X_train.shape}")
    print(f"Target label shape for training: {y_train.shape}\n")
    
    # 2. Create a single dummy test sample
    # This represents a new pair (x1, x2) to be classified
    X_test = np.random.rand(1, 2 * d)
    print(f"Input feature vector shape for a single test sample: {X_test.shape}\n")
    
    # 3. Instantiate, train, and predict with each model
    
    # A. Logistic Regression
    lr_model = LogisticRegression(max_iter=1000)
    lr_model.fit(X_train, y_train)
    lr_prediction = lr_model.predict(X_test)[0]
    print(f"Logistic Regression is SUITABLE. It trained and made a prediction: {lr_prediction}")

    # B. Random Forest
    rf_model = RandomForestClassifier()
    rf_model.fit(X_train, y_train)
    rf_prediction = rf_model.predict(X_test)[0]
    print(f"Random Forest is SUITABLE. It trained and made a prediction: {rf_prediction}")

    # C. Cosine KNN
    # We specify metric='cosine' to use cosine distance
    knn_model = KNeighborsClassifier(metric='cosine')
    knn_model.fit(X_train, y_train)
    knn_prediction = knn_model.predict(X_test)[0]
    print(f"Cosine KNN is SUITABLE. It trained and made a prediction: {knn_prediction}")
    
    print("\n---")
    print("Conclusion: As demonstrated, all listed models are standard classifiers capable of learning from concatenated embedding vectors. Therefore, none of them are unsuitable.")

if __name__ == '__main__':
    demonstrate_model_suitability()