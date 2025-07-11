import numpy as np
from sklearn.linear_model import Lasso

# Plan:
# 1. The user's question is about the strict equivalence between two forms of Lasso:
#    - Penalized (Lagrangian) form: minimize RSS + lambda * L1_norm(beta)
#    - Constrained form: minimize RSS subject to L1_norm(beta) <= t
# 2. This equivalence generally holds in convex optimization but can fail if the solution
#    to the penalized problem is not unique.
# 3. Non-uniqueness can occur, for instance, when the number of predictors (p) is greater
#    than the number of samples (n). In this case, it's possible for different solutions
#    for the same `lambda` to exist and have different L1 norms.
# 4. If solutions for the same `lambda` have different L1 norms (`t`), they cannot all be
#    solutions to the same constrained problem, breaking the strict equivalence.
# 5. Therefore, the statement is "False in general," but it holds in many well-behaved
#    practical cases, making it "often true in practice." This corresponds to option E.
# 6. The following code will illustrate the relationship by solving the penalized form for
#    a specific `lambda` and showing the corresponding constraint `t`, demonstrating the
#    numbers involved in the claimed equivalence.

# This code demonstrates the relationship between the penalty `lambda` (alpha in scikit-learn)
# and the constraint `t` (the L1 norm of the coefficients) in Lasso regression.

# 1. Generate synthetic data with p > n, a case where solution non-uniqueness can occur.
np.random.seed(42)
n_samples, n_features = 50, 100
X = np.random.randn(n_samples, n_features)

# Create a sparse true coefficient vector
n_true_coeffs = 10
true_coef = np.zeros(n_features)
true_coef_indices = np.random.choice(n_features, n_true_coeffs, replace=False)
true_coef[true_coef_indices] = np.random.choice([-5, 5], size=n_true_coeffs)

# Generate target variable y from X, true_coef, and some noise
y = X @ true_coef + np.random.randn(n_samples) * 0.5

# Center the data to handle the intercept implicitly (it becomes zero).
X -= np.mean(X, axis=0)
y -= np.mean(y)

# 2. Solve the penalized Lasso problem for a chosen value of lambda.
# `(α̂, 𝛽̂) = argmin ∑ᵢ(yᵢ — α — ∑ⱼβⱼxᵢⱼ)² + λ∑ⱼ |𝛽ⱼ|`
# We pick a specific lambda (called `alpha` in scikit-learn).
chosen_lambda = 0.5

lasso_penalized = Lasso(alpha=chosen_lambda, max_iter=10000, tol=1e-6)
lasso_penalized.fit(X, y)

# This is the solution `𝛽̂` from the penalized problem
beta_hat = lasso_penalized.coef_

# 3. From this solution, find the corresponding constraint value `t`.
# `t = ∑ⱼ |𝛽̂ⱼ|`
t_constraint = np.sum(np.abs(beta_hat))

# 4. Display the results, illustrating the claimed equivalence with the computed numbers.
# The equivalence states that `beta_hat` should also solve the constrained problem:
# `(α̂, 𝛽̂) = argmin ∑ᵢ(yᵢ — α — ∑ⱼβⱼxᵢⱼ)² subject to ∑ⱼ |𝛽ⱼ| ≤ t`

print("Demonstrating the relationship between the two Lasso formulations.")
print(f"\n1. Solved the Penalized Form with lambda = {chosen_lambda}")
print("   - This minimizes: RSS + lambda * (L1 norm of beta)")
print(f"   - The resulting coefficient vector 'beta_hat' has an L1 norm of: {t_constraint:.4f}\n")


print("2. The corresponding Constrained Form problem is:")
print("   - Minimize: RSS")
print(f"   - Subject to: (L1 norm of beta) <= {t_constraint:.4f}\n")


print("The core of the question is whether the solution set is ALWAYS identical for some (lambda, t) pair.")
print("The answer is no, because when the penalized problem has multiple solutions, those solutions")
print("can have different L1 norms, so they can't all solve the same constrained problem.\n")


# Per the prompt, output the numbers in the final equation.
print("--- Final Equation Numbers ---")
y_pred = X @ beta_hat
rss = np.sum((y - y_pred)**2)
l1_penalty_value = chosen_lambda * t_constraint

print(f"For the PENALIZED problem with lambda = {chosen_lambda}:")
print(f"  argmin( RSS + {chosen_lambda} * sum(|beta_j|) )")
print(f"  At the solution, RSS = {rss:.4f} and the penalty term = {l1_penalty_value:.4f}")
print(f"  Total objective value = {rss + l1_penalty_value:.4f}")

print("\nFor the CONSTRAINED problem:")
print(f"  argmin( RSS )  subject to sum(|beta_j|) <= {t_constraint:.4f}")
print(f"  The equivalence implies the minimal RSS under this constraint is {rss:.4f}")
