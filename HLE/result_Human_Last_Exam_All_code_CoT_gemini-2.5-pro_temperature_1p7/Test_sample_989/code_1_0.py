import numpy as np
from sklearn.svm import SVC
from sklearn.datasets import make_blobs

# This script provides a practical demonstration for statement D:
# "It is possible to add or move data points and not affect the decision boundary at all,
# as long as they're interior points"
# We will create an SVM, then add a point far from the boundary and show that
# the model's decision function does not change.

# 1. Create an initial dataset and train an SVM
X, y = make_blobs(n_samples=50, centers=2, random_state=4, cluster_std=1.0)
svm_initial = SVC(kernel='linear', C=100)
svm_initial.fit(X, y)

# 2. Store the coefficients of the initial model's decision boundary
initial_w = svm_initial.coef_[0]
initial_b = svm_initial.intercept_[0]

# 3. Add a new data point that is clearly an "interior point"
# This point is correctly classified and far from the separating hyperplane.
new_point = np.array([[-10, -10]])
new_label = np.array([0]) # Assuming the cluster around (-10,-10) is class 0

X_updated = np.vstack([X, new_point])
y_updated = np.hstack([y, new_label])

# 4. Retrain the SVM on the dataset with the new point
svm_updated = SVC(kernel='linear', C=100)
svm_updated.fit(X_updated, y_updated)

# 5. Get the coefficients of the updated model
updated_w = svm_updated.coef_[0]
updated_b = svm_updated.intercept_[0]

# 6. Print and compare the results
print("--- Initial Model ---")
print("The decision function is defined by w * x + b = 0")
print(f"w1 = {initial_w[0]}")
print(f"w2 = {initial_w[1]}")
print(f"b  = {initial_b}")
print(f"Equation: {initial_w[0]:.4f} * x1 + {initial_w[1]:.4f} * x2 + {initial_b:.4f} = 0")
print(f"Number of support vectors: {len(svm_initial.support_vectors_)}")

print("\n--- Model After Adding Interior Point ---")
print("A new point at (-10, -10) was added, far from the decision boundary.")
print(f"w1 = {updated_w[0]}")
print(f"w2 = {updated_w[1]}")
print(f"b  = {updated_b}")
print(f"Equation: {updated_w[0]:.4f} * x1 + {updated_w[1]:.4f} * x2 + {updated_b:.4f} = 0")
print(f"Number of support vectors: {len(svm_updated.support_vectors_)}")

# Final check
if np.allclose(initial_w, updated_w) and np.allclose(initial_b, updated_b):
    print("\nConclusion: The decision boundary is unchanged, demonstrating that interior points do not affect the SVM solution.")
else:
    print("\nConclusion: The decision boundary changed.")
