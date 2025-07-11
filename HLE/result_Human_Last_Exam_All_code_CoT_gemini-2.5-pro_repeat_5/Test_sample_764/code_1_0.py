import numpy as np
from scipy.optimize import minimize

# We will investigate the equivalence of the two Lasso formulations with a simple
# counterexample where the number of predictors p=2 is greater than the number
# of observations n=1.
# Let y = [1] and X = [[1, 1]]. We will ignore the intercept for simplicity.
# The residual sum of squares (RSS) is (y - (b1*x1 + b2*x2))^2,
# which simplifies to (1 - b1 - b2)^2.

# --- Problem 1: Constrained Lasso ---
# argmin (1 - b1 - b2)^2  subject to |b1| + |b2| <= t

# Let's choose t = 1. The constraint is |b1| + |b2| <= 1.
# The RSS is minimized when its value is 0, which occurs when b1 + b2 = 1.
# We need to find the points that satisfy both b1 + b2 = 1 and |b1| + |b2| <= 1.
# If b1 + b2 = 1 and we assume b1>=0, b2>=0, then |b1|+|b2| = b1+b2 = 1.
# The constraint is met exactly on the boundary.
# So, the solution set for the constrained problem with t=1 is the line segment
# connecting (1, 0) and (0, 1).
# S1(t=1) = {(b1, b2) | b1 + b2 = 1, b1 >= 0, b2 >= 0}.
# This is an infinite set of solutions.

print("--- Analysis of the Constrained Problem (t=1) ---")
print("Objective function: (1 - b1 - b2)^2")
print("Constraint: |b1| + |b2| <= 1")
print("The minimum objective value is 0.")
print("This minimum is achieved for any point (b1, b2) on the line b1 + b2 = 1.")
print("The intersection of the solution space (b1+b2=1) and the constraint region (|b1|+|b2|<=1) is the line segment from (1, 0) to (0, 1).")
print("Thus, the solution set S1(t=1) is this entire line segment, containing infinite solutions.\n")

# Let's verify one point on this segment, e.g., (0.5, 0.5)
b_constrained_example = np.array([0.5, 0.5])
rss_val = (1 - b_constrained_example.sum())**2
l1_norm = np.sum(np.abs(b_constrained_example))
print(f"Example solution for t=1: beta = {b_constrained_example}")
print(f"RSS at this point: {rss_val:.4f}")
print(f"L1 norm at this point: {l1_norm:.4f}, which satisfies the constraint.\n")


# --- Problem 2: Penalized Lasso ---
# argmin (1 - b1 - b2)^2 + lambda * (|b1| + |b2|)

# Let's analyze this problem. For any lambda > 0.
# Let u = b1 + b2. We are minimizing (1 - u)^2 + lambda * (|b1| + |b2|).
# For a fixed u = b1 + b2, the L1 norm |b1| + |b2| is minimized when one
# coefficient is u and the other is 0. This is a key property of the L1 penalty.
# So, the solution will have the form (u, 0) or (0, u).
# The problem reduces to minimizing f(u) = (1 - u)^2 + lambda * |u|.
# By solving df/du = 0, we find the solution is u = max(0, 1 - lambda/2).
# This means for any lambda in (0, 2), the solution u is in (0, 1).
# The solution set S2(lambda) will be {(u, 0), (0, u)}. This is a set of TWO points.

# For example, if we seek a solution matching the constrained case where t=1, we might
# try to find a lambda that gives a solution with L1 norm of 1.
# This happens when u=1, which corresponds to lambda -> 0.
# The solution set for lambda -> 0+ is {(1, 0), (0, 1)}.
# This two-point set is NOT the same as S1(t=1), which was the entire line segment.

print("--- Analysis of the Penalized Problem ---")
print("Objective function: (1 - b1 - b2)^2 + lambda * (|b1| + |b2|)")
print("For any lambda > 0, the solution set has at most two points: (u, 0) and (0, u).")

def penalized_loss(b, lam):
    return (1 - b.sum())**2 + lam * np.sum(np.abs(b))

# Let's use lambda = 1. Analytically, u = 1 - 1/2 = 0.5.
# The solution set S2(lambda=1) is {(0.5, 0), (0, 0.5)}.
res = minimize(penalized_loss, x0=[0.1, 0.2], args=(1.0,), tol=1e-9)
b_penalized = res.x
print(f"\nNumerical solution for lambda=1: beta = [{b_penalized[0]:.6f}, {b_penalized[1]:.6f}]")
# Note: The optimizer finds one of the two solutions.

sol_from_s1 = np.array([0.5, 0.5]) # A point from the constrained solution set
loss_s1_point = penalized_loss(sol_from_s1, 1.0)
loss_s2_point = penalized_loss(b_penalized, 1.0)
print(f"Loss at the penalized solution ({b_penalized[0]:.1f}, {b_penalized[1]:.1f}): {loss_s2_point:.4f}")
print(f"Loss at the point (0.5, 0.5) from S1(t=1): {loss_s1_point:.4f}")
print("The minimum for the penalized problem is at a corner, not the center of the segment.\n")

print("--- Conclusion ---")
print("The solution set for the constrained problem with t=1 is a line segment (infinite solutions).")
print("The solution set for the penalized problem is always a finite set of points (for lambda > 0).")
print("Since the solution sets are fundamentally different, the two formulations are not strictly equivalent.")
print("Therefore, the original statement is false.")
print("Final equation for one point in the penalized solution set with lambda=1: b1 + b2 = 0.5")
print("Example solution values: b1 = 0.5, b2 = 0.0")
