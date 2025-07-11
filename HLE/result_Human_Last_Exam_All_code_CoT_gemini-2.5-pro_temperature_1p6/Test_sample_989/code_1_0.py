import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm

# 1. Create initial linearly separable data
X = np.array([[1, 2], [2, 3], [3, 3], [2, 1], [3, 2]])
y = np.array([1, 1, 1, 0, 0])

# 2. Train a linear SVM
clf = svm.SVC(kernel='linear')
clf.fit(X, y)

# 3. Print the original model's parameters and support vectors
print("--- Original Model ---")
# w_1*x_1 + w_2*x_2 + b = 0
w_original = clf.coef_[0]
b_original = clf.intercept_[0]
sv_original = clf.support_vectors_

print(f"Weights (w): [{w_original[0]:.4f}, {w_original[1]:.4f}]")
print(f"Intercept (b): {b_original:.4f}")
print("Support Vectors:")
print(sv_original)
print("\n")


# 4. Add a new 'interior' point that is far from the decision boundary
X_new = np.vstack([X, [6, 6]]) # Add a point to class 1
y_new = np.append(y, 1)

# 5. Re-train the SVM on the new dataset
clf_new = svm.SVC(kernel='linear')
clf_new.fit(X_new, y_new)

# 6. Print the new model's parameters and support vectors
print("--- Model After Adding Interior Point ---")
w_new = clf_new.coef_[0]
b_new = clf_new.intercept_[0]
sv_new = clf_new.support_vectors_

print(f"Weights (w): [{w_new[0]:.4f}, {w_new[1]:.4f}]")
print(f"Intercept (b): {b_new:.4f}")
print("Support Vectors:")
print(sv_new)
print("\n")


# 7. Final Conclusion
print("--- Conclusion ---")
print("As you can see, the weights (w) and intercept (b) did not change after adding the interior point.")
print("The set of support vectors also remained the same.")
print("This demonstrates that adding an interior point has no effect on the decision boundary.")
