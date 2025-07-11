import numpy as np
from sklearn.linear_model import Lasso
from sklearn.linear_model import LinearRegression

# Plan: Demonstrate that the mapping between LASSO's parameters t and λ is not always 1-to-1.

# 1. Generate some synthetic data for a regression problem.
np.random.seed(42)
n_samples, n_features = 50, 20
X = np.random.randn(n_samples, n_features)
# Create a sparse ground truth for coefficients
true_coef = 3 * np.random.randn(n_features)
inds = np.arange(n_features)
np.random.shuffle(inds)
true_coef[inds[5:]] = 0  # Make most coefficients zero
y = X @ true_coef + np.random.randn(n_samples) * 0.5

# For simplicity, we center the data. sklearn's Lasso handles the intercept automatically.
X_centered = X - np.mean(X, axis=0)
y_centered = y - np.mean(y)

# 2. Solve the penalized LASSO problem for a wide range of λ (alpha in sklearn).
# We go from a small λ (close to OLS) to a large λ (all coefficients zero).
alphas = np.logspace(-3, 2, 300)
l1_norms = []

for a in alphas:
    # Solve argmin ||y - Xβ||² + λ ||β||₁
    lasso = Lasso(alpha=a, fit_intercept=False, tol=1e-8, max_iter=10000)
    lasso.fit(X_centered, y_centered)
    # Calculate t = ||β||₁ for the solution
    l1_norms.append(np.sum(np.abs(lasso.coef_)))

# 3. Analyze the mapping at the boundaries.

# Case 1: Large t / λ = 0 (Ordinary Least Squares)
# The constrained problem for any t >= ||β_OLS||₁ gives the OLS solution.
# This corresponds to λ=0 in the penalized version.
ols = LinearRegression(fit_intercept=False)
ols.fit(X_centered, y_centered)
t_ols = np.sum(np.abs(ols.coef_))

print("Demonstration of the relationship between parameters t and λ in LASSO")
print("=" * 60)
print("CASE 1: Many-to-one mapping (t -> λ)")
print(f"The Ordinary Least Squares (OLS) solution corresponds to λ = 0.")
print(f"The L1 norm (our 't') for the OLS solution is: {t_ols:.4f}")
print("In the constrained formulation, any t >= {t_ols:.4f} results in this same OLS solution.")
print("For example, a constraint of t = {t_ols + 1:.4f} and t = {t_ols + 10:.4f} both map to λ = 0.")
print("This is a many-to-one mapping from t to λ.\n")


# Case 2: t = 0 / Large λ
# The constrained problem for t=0 forces all coefficients to be zero.
# This corresponds to any sufficiently large λ in the penalized version.
# Find the first alpha where the L1 norm is effectively zero.
is_zero = np.isclose(l1_norms, 0)
if np.any(is_zero):
    lambda_min_for_zero_solution = alphas[is_zero].min()
else:
    lambda_min_for_zero_solution = alphas[-1] # Should not happen with large alpha range

print("CASE 2: One-to-many mapping (t -> λ)")
print("A constraint of t = 0 forces all coefficients to be zero (β=0).")
print(f"In the penalized form, this β=0 solution is achieved for any λ >= {lambda_min_for_zero_solution:.4f}.")
print(f"For λ = {lambda_min_for_zero_solution:.4f}, the L1 norm t is {l1_norms[list(alphas).index(lambda_min_for_zero_solution)]:.6f}")
large_lambda = alphas[-1]
print(f"For λ = {large_lambda:.4f}, the L1 norm t is {l1_norms[-1]:.6f}")
print("This is a one-to-many mapping from t=0 to λ.\n")


# 4. Conclusion
print("CONCLUSION:")
print("While the two LASSO formulations produce the same set of possible solutions (so the equivalence is 'true in general'),")
print("the relationship between their parameters, t and λ, is not a strict one-to-one mapping.")
print("Therefore, stating the equivalence is 'always strictly true' is an overstatement.")
print("=" * 60)
