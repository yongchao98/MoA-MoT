import numpy as np
from sklearn.linear_model import Lasso
import matplotlib.pyplot as plt

# Plan:
# 1. Generate synthetic data where the number of features (p) is greater than the number of samples (n).
# 2. Define a range of regularization parameters, lambda (alpha in scikit-learn).
# 3. For each lambda, fit a LASSO model and calculate the L1 norm of the resulting coefficients (this is 't').
# 4. Plot 't' vs. lambda to observe the relationship.
# 5. Identify and comment on any "plateau" regions, where 't' is nearly constant for a range of lambda values.
#    This plateau demonstrates the breakdown of a strict one-to-one mapping.

# Step 1: Generate synthetic data with p > n
np.random.seed(42)
n_samples, n_features = 50, 100 # p > n case
X = np.random.randn(n_samples, n_features)

# Generate a sparse true coefficient vector to create a ground truth
true_coef = np.zeros(n_features)
true_coef[10:20] = 5
y = X @ true_coef + np.random.randn(n_samples) * 0.5

# Step 2: Define a range of lambda (alpha) values
# We use a logarithmic scale to cover a wide range of penalties.
alphas = np.logspace(-4, 2, 300)

# Step 3: Fit LASSO for each alpha and compute the L1 norm 't'
l1_norms = []
for a in alphas:
    # We use scikit-learn's Lasso, where alpha is the regularization parameter lambda.
    # The objective is slightly different (scaled), but the solution path behavior is the same.
    lasso = Lasso(alpha=a, fit_intercept=True, max_iter=10000, tol=1e-5)
    lasso.fit(X, y)
    t = np.sum(np.abs(lasso.coef_))
    l1_norms.append(t)

# --- Output and Explanation ---

print("Demonstrating the relationship between the LASSO constraint 't' and penalty 'lambda'.")
print(f"We simulated a case with n={n_samples} samples and p={n_features} features (where p > n).\n")

# Find the approximate start of the plateau for illustrative purposes.
# The plateau forms for small alpha, where the L1 norm becomes large.
# We'll look at the last part of our collected data (corresponding to small alphas).
end_l1_norm = l1_norms[-1]
plateau_start_index = np.where(np.array(l1_norms) > end_l1_norm * 0.99)[0][0]
plateau_alpha = alphas[plateau_start_index]
plateau_t = l1_norms[plateau_start_index]

print(f"In our simulation, for lambda (alpha) values smaller than ~{plateau_alpha:.4f}, the L1 norm 't' becomes nearly constant.")
print("Example values from the plateau region:")
# We reverse the lists to show lambda decreasing and 't' increasing towards the plateau.
for i in range(5):
    idx = len(alphas) - 1 - i * 5
    alpha_val = alphas[idx]
    t_val = l1_norms[idx]
    print(f"  For lambda = {alpha_val:8.4f}, the corresponding t (L1 norm) is {t_val:8.4f}")

print("\nAs you can see, a range of different 'lambda' values all produce almost the same 't' value.")
print("This shows that there isn't always a unique lambda for a given t.")
print("Therefore, the equivalence is not 'always strictly true'.")

# Step 4: Plot t vs. lambda
plt.figure(figsize=(12, 7))
plt.plot(alphas, l1_norms, marker='.')
plt.xscale('log')
plt.xlabel('Regularization Strength λ (alpha)')
plt.ylabel('L1 Norm of Coefficients t = Σ|β|')
plt.title('t vs. λ for LASSO with More Features than Samples (p > n)')
plt.grid(True)
# Invert x-axis so lambda=0 is on the left, making the path easier to interpret
plt.gca().invert_xaxis()
# Highlight the plateau region
plt.axvspan(alphas[-1], plateau_alpha, color='red', alpha=0.2, label='Plateau Region (Non-unique λ for a given t)')
plt.legend()
plt.show()
