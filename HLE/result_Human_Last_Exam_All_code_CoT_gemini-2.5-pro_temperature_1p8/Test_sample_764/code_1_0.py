import numpy as np
from sklearn.linear_model import LinearRegression

# Plan:
# 1. Define a simple regression problem with p > n and collinearity,
#    where the Ordinary Least Squares (OLS) solution is not unique.
# 2. Analyze the solution set for the penalized LASSO formulation with lambda = 0 (which is OLS).
# 3. Analyze the solution set for the constrained LASSO formulation with t chosen to be the
#    minimum L1 norm among all OLS solutions.
# 4. Show that these two solution sets are not identical, which proves the equivalence is not strictly true.

# Step 1: Define the problem
# Let's use a simple problem: p=2 predictors, n=1 observation.
# The predictors are perfectly collinear.
X = np.array([[1, 1]])
y = np.array([1])

print("Demonstrating the non-equivalence of LASSO formulations.")
print("="*60)
print(f"Problem setup: X = {X}, y = {y}")
print("The model is y = β₁*x₁ + β₂*x₂. Here, x₁=1, x₂=1.")
print("The Residual Sum of Squares (RSS) is (y - (β₁*x₁ + β₂*x₂))² = (1 - (β₁ + β₂))²\n")

# Step 2: Analyze the penalized formulation with lambda = 0 (OLS)
print("--- Penalized Formulation (with λ = 0) ---")
print("This is equivalent to Ordinary Least Squares (OLS).")
print("Objective: minimize (1 - (β₁ + β₂))²")
print("The minimum RSS is 0, which is achieved by any pair of coefficients (β₁, β₂) as long as their sum is 1.")
print("The set of solutions, let's call it B(λ=0), is the entire line defined by the equation:")
print("β₁ + β₂ = 1")
print("For example, (1, 0), (0, 1), (2, -1), (0.5, 0.5) are all valid OLS solutions.\n")

# Use scikit-learn to find one such solution.
# Note: Lasso(alpha=0) gives a warning and defaults to LinearRegression.
ols_model = LinearRegression(fit_intercept=False)
ols_model.fit(X, y)
b_ols = ols_model.coef_
print(f"Scikit-learn's LinearRegression finds one solution from this infinite set: β = {np.round(b_ols, 2)}")
print("This specific solution is an artifact of the algorithm, but the true solution set is the infinite line.\n")


# Step 3: Analyze the constrained formulation
print("--- Constrained Formulation ---")
print("Objective: minimize (1 - (β₁ + β₂))² subject to |β₁| + |β₂| ≤ t")
print("First, we must find the OLS solutions that have the minimum L1 norm (|β₁|+|β₂|).")
print("We need to minimize |β₁| + |β₂| under the condition that b1 + b2 = 1.")
print("Let β₂ = 1 - β₁. We minimize f(β₁) = |β₁| + |1 - β₁|.")
print("This function f(β₁) is equal to 1 for any β₁ in [0, 1], and greater than 1 otherwise.")
print("The minimum L1 norm is 1. We therefore set our constraint parameter t = 1.\n")

print(f"Now, let's solve the constrained problem for t = 1:")
print("To get the minimum RSS of 0, we still need b1 + b2 = 1.")
print("The solution must satisfy BOTH b₁ + b₂ = 1 AND |β₁| + |β₂| ≤ 1.")
print("The set of (β₁, β₂) where β₁ + β₂ = 1 and |β₁| + |β₂| ≤ 1 is the line segment from (1,0) to (0,1).")
print("So the solution set, A(t=1), is defined by the equations:")
print("β₁ + β₂ = 1, where 0 ≤ β₁ ≤ 1\n")

# Step 4: Compare the solution sets
print("--- Conclusion ---")
print("The solution set for the penalized problem with λ=0 is an infinite line.")
print("The solution set for the constrained problem with the corresponding t=1 is a finite line segment.")
print("\nClearly, A(t=1) is a proper subset of B(λ=0). The sets are NOT identical.")
print("Therefore, the strict equivalence is false in general, though it often holds in practice under ideal conditions (e.g., unique solutions).")