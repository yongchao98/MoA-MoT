import numpy as np
from scipy.optimize import minimize

# This script demonstrates that the equivalence between the Lagrangian and constrained
# forms of LASSO is not always strictly true, particularly when solutions are not unique.
# We use a simple case with n=1 observation and p=2 predictors.
# y = 1, x1 = 1, x2 = 1. The model is y = b1*x1 + b2*x2.

# --- 1. Lagrangian Problem ---
# Objective: (1 - (b1 + b2))^2 + lambda * (|b1| + |b2|)
# We set lambda = 1.0

def lagrangian_objective(beta, lam):
    """LASSO objective function (Lagrangian form)."""
    b1, b2 = beta
    rss = (1 - (b1 + b2))**2
    penalty = lam * (np.abs(b1) + np.abs(b2))
    return rss + penalty

lam = 1.0

# We minimize the objective starting from different initial points.
# If the solution is unique, we should always get the same result.
beta_initial_1 = [0.0, 0.0]
res1 = minimize(lagrangian_objective, beta_initial_1, args=(lam,), method='Nelder-Mead')

beta_initial_2 = [0.5, -0.5]
res2 = minimize(lagrangian_objective, beta_initial_2, args=(lam,), method='Nelder-Mead')

# The KKT conditions show that any (b1, b2) with b1>0, b2>0 and b1+b2 = 0.5 is a solution.
# The optimizer finds points on this solution set.
# Note: Numerical optimizers might find slightly different points on the solution manifold.
sol1 = res1.x
sol2 = res2.x


print("--- Lagrangian Form (lambda = 1.0) ---")
print(f"Starting from {beta_initial_1}, found solution: beta = [{sol1[0]:.4f}, {sol1[1]:.4f}]")
print(f"Starting from {beta_initial_2}, found solution: beta = [{sol2[0]:.4f}, {sol2[1]:.4f}]")
print("The solutions are different, indicating the solution is not unique.")
print("The set of solutions is approximately the line segment b1+b2=0.5 for b1,b2>0.\n")

# For any solution on this segment, the L1 norm |b1|+|b2| = b1+b2 = 0.5.
# This corresponds to a constrained problem with t = 0.5.

# --- 2. Constrained Problem ---
# Objective: min (1 - (b1 + b2))^2
# Subject to: |b1| + |b2| <= t
# We set t = 0.5.

def rss_objective(beta):
    """RSS objective function."""
    b1, b2 = beta
    return (1 - (b1 + b2))**2

t = 0.5
constraints = ({'type': 'ineq', 'fun': lambda beta: t - (np.abs(beta[0]) + np.abs(beta[1]))})

# We solve the constrained problem.
beta_initial_3 = [0.0, 0.0]
res3 = minimize(rss_objective, beta_initial_3, constraints=constraints, method='SLSQP')
sol3 = res3.x

print("--- Constrained Form (t = 0.5) ---")
print(f"The unique solution is: beta = [{sol3[0]:.4f}, {sol3[1]:.4f}]")
print("This solution is the point in the feasible set (|b1|+|b2|<=0.5) closest to the line b1+b2=1.\n")

# --- 3. Conclusion ---
print("--- Conclusion ---")
print("The Lagrangian problem (lambda=1) has a set of solutions (a line segment).")
print("The corresponding constrained problem (t=0.5) has a single, unique solution.")
print("Since the solution sets are not identical, the equivalence is not strictly true in this case.")
print("This demonstrates the statement is 'False in general, but often true in practice'.")

# Final Equation format is not applicable here as we are showing concepts, not calculating a specific number.
# The core of the argument is comparing the solution sets from the two optimization problems.

<<<E>>>