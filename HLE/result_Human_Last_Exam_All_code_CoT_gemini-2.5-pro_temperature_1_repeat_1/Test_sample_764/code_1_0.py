import numpy as np
from sklearn.linear_model import LinearRegression

# 1. Explain the overall logic
print("We will investigate the equivalence between two forms of Lasso regression.")
print("Form 1 (Constrained): min RSS(β) subject to ∑|βⱼ| ≤ t")
print("Form 2 (Penalized):   min RSS(β) + λ∑|βⱼ|\n")
print("The question is whether for any `t`, there's a `λ` (and vice-versa) that yields the same solution.")
print("We will show this is generally true, but not strictly true, by demonstrating a case where the equivalence breaks down.\n")
print("Our plan is to show that for a budget `t` larger than the L1 norm of the OLS solution,")
print("the equivalence does not strictly hold.\n")

# 2. Setup a simple regression problem
# Using np.newaxis to make them column vectors for sklearn
X = np.array([1, 2, 3, 4, 5]).reshape(-1, 1)
y = np.array([2.5, 3.5, 4.0, 5.5, 6.0])

# 3. Solve for OLS solution and its L1 norm (t_OLS)
print("--- Step 1: Find the OLS Solution ---")
# Using a simplified model with no intercept for clarity
ols_model = LinearRegression(fit_intercept=False)
ols_model.fit(X, y)
beta_ols = ols_model.coef_
t_ols = np.sum(np.abs(beta_ols))

print(f"The OLS solution is β_OLS = {beta_ols[0]:.4f}")
print(f"The L1 norm of the OLS solution is t_OLS = ∑|β_OLS| = {t_ols:.4f}\n")

# 4. Analyze the constrained problem for a case where t > t_OLS
print("--- Step 2: Analyze the Constrained Problem ---")
t_user = t_ols + 2.0  # Choose a t larger than t_OLS
print(f"Let's consider the constrained problem with a budget t = {t_user:.4f}")
print(f"Since our chosen t ({t_user:.4f}) is greater than t_OLS ({t_ols:.4f}), the constraint ∑|β| ≤ {t_user:.4f} is inactive.")
print("The solution that minimizes the RSS is the unconstrained OLS solution itself.")
solution_constrained = beta_ols
print(f"Therefore, the solution to the constrained problem is β_constrained = {solution_constrained[0]:.4f}\n")

# 5. Analyze the corresponding penalized problem
print("--- Step 3: Analyze the Penalized Problem and Check Equivalence ---")
print(f"Now, we look for a penalized problem 'min RSS(β) + λ∑|β|' that gives the solution β = {solution_constrained[0]:.4f}.")
print("The only way to get the unconstrained OLS solution from the penalized form is to have no penalty.")
lambda_equiv = 0.0
print(f"This requires the penalty term λ to be {lambda_equiv:.1f}.\n")

print("For the penalized problem with λ = 0.0, the solution is the OLS solution.")
solution_penalized = beta_ols # The solution for lambda=0 is beta_ols by definition
t_from_penalized = np.sum(np.abs(solution_penalized))

print(f"The solution β_penalized for λ=0 is indeed {solution_penalized[0]:.4f}.")
print(f"The L1 norm (the 't' value) corresponding to this solution is ∑|β_penalized| = {t_from_penalized:.4f}.\n")


# 6. Conclusion
print("--- Conclusion: Comparing the 't' values ---")
print("We started with the constrained problem using a specific budget `t_user`.")
print("We found the solution `β_constrained`.")
print("We then found the `λ` that produces this same solution in the penalized problem.")
print("Finally, we compare the original budget `t_user` with the L1 norm derived from the equivalent penalized problem.")

# Outputting the numbers for the final comparison
print("\nFinal Comparison:")
print(f"Original budget from constrained problem: t_user           = {t_user:.4f}")
print(f"L1 norm from the equivalent penalized problem: t_from_penalized = {t_from_penalized:.4f}")

if not np.isclose(t_user, t_from_penalized):
    print("\nResult: The two 't' values are NOT equal.")
    print("This demonstrates that while we can find a λ that gives the right β, the 't' values do not match.")
    print("Therefore, the equivalence between the two forms is not strictly true for all possible values of t.")
else:
    print("\nResult: The two 't' values are equal, which contradicts the expected outcome for this edge case.")

print("\nThis illustrates the point in option K: The equivalence is 'True in general, but not always strictly true'.")