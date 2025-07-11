import numpy as np

# We investigate the strict equivalence of two Lasso formulations:
# 1. Constrained: argmin RSS(β) subject to ||β||_1 <= t
# 2. Penalized:   argmin RSS(β) + λ * ||β||_1
#
# We will show that this equivalence is NOT strictly true by providing a counterexample where p > n.

print("--- Step 1: Define a Counterexample (p > n) ---")
# Let's consider a simple regression problem with n=1 observation and p=2 predictors.
# We'll set the intercept α=0 for simplicity.
# y = [1]
# X = [[1, 0]]
# The regression model is y_i = β₁*xᵢ₁ + β₂*xᵢ₂
# For our single observation (i=1): 1 = β₁*1 + β₂*0, which means 1 = β₁.
# The Residual Sum of Squares (RSS) is (y - Xβ)² = (1 - β₁)²
print("Consider a problem with y=[1] and X=[[1, 0]].")
print(f"The RSS to minimize is (1 - \u03B2\u2081)\u00B2.")
print("Notice that RSS depends only on \u03B2\u2081 and is minimized (to 0) when \u03B2\u2081 = 1.")
print("\u03B2\u2082 can be any value without changing the RSS, so the OLS solution is not unique.")
print("-" * 50)


print("--- Step 2: Analyze the Constrained Problem ---")
# Problem: argmin (1 - β₁)² subject to |β₁| + |β₂| <= t
print(f"Let's solve the constrained problem for a budget t = 1.5.")
print(f"We want to minimize (1 - \u03B2\u2081)\u00B2 subject to |\u03B2\u2081| + |\u03B2\u2082| \u2264 1.5.")
print("To achieve the minimum RSS of 0, we must have \u03B2\u2081 = 1.")
print("With \u03B2\u2081 = 1, the constraint becomes |1| + |\u03B2\u2082| \u2264 1.5, which means |\u03B2\u2082| \u2264 0.5.")
print("So, the set of solutions for the constrained problem is: {(1, \u03B2\u2082) | -0.5 \u2264 \u03B2\u2082 \u2264 0.5}.")
print(f"Let's pick ONE solution from this set: (\u03B2\u0302\u2081, \u03B2\u0302\u2082) = (1, 0.5).")
print("-" * 50)


print("--- Step 3: Analyze the Penalized Problem ---")
# Problem: argmin (1 - β₁)² + λ * (|β₁| + |β₂|)
print("Now let's analyze the penalized problem for any λ > 0.")
print(f"Objective: Minimize (1 - \u03B2\u2081)\u00B2 + \u03BB(|\u03B2\u2081| + |\u03B2\u2082|).")
print("To minimize this, we must minimize each part with respect to the variables.")
print(f"The term \u03BB*|\u03B2\u2082| is non-negative and is minimized (to 0) only when \u03B2\u2082 = 0.")
print("Therefore, any solution to the penalized problem (for λ > 0) MUST have \u03B2\u2082 = 0.")
print("-" * 50)


print("--- Step 4: The Contradiction ---")
beta1_hat_c, beta2_hat_c = 1, 0.5
print(f"Our chosen solution from the constrained problem was (\u03B2\u0302\u2081, \u03B2\u0302\u2082) = ({beta1_hat_c}, {beta2_hat_c}).")
print("Our analysis of the penalized problem showed that any solution must have \u03B2\u2082 = 0.")
print(f"Since {beta2_hat_c} \u2260 0, our constrained solution cannot be a solution to the penalized problem for any \u03BB > 0.")
print("\nThis counterexample shows that the solution sets of the two problems are not always identical.")
print("The statement of strict equivalence is therefore false.")
print("\nThis failure happens in cases like p > n or with perfect collinearity. In many well-behaved practical cases (p <= n, no collinearity), the equivalence holds.")
print("This makes 'False in general, but often true in practice' the best description.")
