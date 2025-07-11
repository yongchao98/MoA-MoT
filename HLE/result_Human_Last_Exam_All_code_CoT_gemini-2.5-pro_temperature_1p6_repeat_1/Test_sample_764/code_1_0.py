import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import lasso_path
from sklearn.datasets import make_regression

# Step 1 & 2: Generate synthetic data and compute the Lasso path to demonstrate its continuity.
# The Lasso path shows how the regression coefficients change as the regularization parameter lambda (α) varies.
# This visualization will help in evaluating the answer choices.
X, y = make_regression(n_samples=100, n_features=10, n_informative=5, noise=25, random_state=42)

# The function `lasso_path` computes the coefficients along the path.
# We choose a range of alpha values to trace the path from the OLS solution towards zero.
alphas, coefs, _ = lasso_path(X, y, eps=1e-3, n_alphas=200, fit_intercept=False)

# Step 3: Plot the coefficient paths.
# Each line in the plot represents the value of a single coefficient β_j as a function of α.
plt.figure(figsize=(12, 7))
colors = plt.cm.get_cmap('tab10', X.shape[1])

for i in range(X.shape[1]):
    # We plot |β_j| vs. ||β||₁ to show the evolution relative to the budget 't'.
    # This is an alternative and common way to view the LARS-Lasso path.
    l1_norm = np.sum(np.abs(coefs), axis=0)
    plt.plot(l1_norm, coefs[i, :], color=colors(i), label=f'$\\beta_{{{i+1}}}$')

# Final plot formatting
plt.xlabel("L1 Norm of Coefficients ($t = \\sum_j |\\beta_j|$)")
plt.ylabel("Coefficient Values ($\\beta_j$)")
plt.title("Lasso Coefficient Paths")
plt.legend(loc='upper left')
plt.grid(True)
plt.show()

# Step 4: Explain the conclusion from the plot and theory.
print("The plot demonstrates that the Lasso coefficients (β_j) change continuously as their L1 norm (the budget 't') changes.")
print("This visual evidence contradicts the reasoning in choices A, C, and H, which incorrectly claim the path is not continuous.")
print("\nThe two formulations of Lasso are:")
print("1. Constrained: argmin RSS(β) subject to ∑|β_j| ≤ t")
print("2. Penalized: argmin RSS(β) + λ∑|β_j|")
print("\nThese two forms are equivalent under most practical conditions, forming the basis of the LARS algorithm.")
print("However, the equivalence is not *strictly* universal. It can fail if the solution is not unique (e.g., with perfectly collinear predictors).")
print("In such special cases, the set of minimizers for the constrained problem might not be identical to the set for the penalized problem.")
print("Therefore, the statement 'always strictly true' is false. The best description is that it's 'False in general, but often true in practice'.")
