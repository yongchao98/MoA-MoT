import numpy as np

print("We want to find the normal cone T_F^°(x^*) for the given feasible set F at x* = (2, 0, -1)^T.")

# Step 1: Define the problem specifics
x_star = np.array([2, 0, -1])

# The four inequality constraints g_i(x) <= 0 are defined as functions:
def g1(x):
    return (x[0] - 1)**2 + x[1]**2 - 1
def g2(x):
    return (x[0] - 3)**2 + x[1]**2 - 1
def g3(x):
    return x[2] + 1
def g4(x):
    return -x[2] - 2

constraints = [g1, g2, g3, g4]

# Step 2: Identify active constraints at x*
print("\n--- Step 1: Identify Active Constraints ---")
print(f"We evaluate each constraint function g_i(x) at x* = {x_star.tolist()}:")

g_values = [g(x_star) for g in constraints]
active_indices = []
for i, val in enumerate(g_values):
    status = "Active" if np.isclose(val, 0) else "Inactive"
    if status == "Active":
        active_indices.append(i)
    print(f"g_{i+1}(x*) = {val:.4f}  => Constraint {i+1} is {status}")

print(f"\nThe set of active constraint indices is I(x*) = {{ {', '.join(str(i+1) for i in active_indices)} }}")

# Step 3: Compute gradients of active constraints
print("\n--- Step 2: Compute Gradients of Active Constraints ---")
print("Next, we compute the gradients of the active constraints at x*.")

# Gradients of the constraint functions
def grad_g1(x):
    return np.array([2 * (x[0] - 1), 2 * x[1], 0])
def grad_g2(x):
    return np.array([2 * (x[0] - 3), 2 * x[1], 0])
def grad_g3(x):
    return np.array([0, 0, 1])
# g4 is not active, so its gradient is not needed for the cone.
active_grad_funcs = [grad_g1, grad_g2, grad_g3]
grad_vectors = []
for i, grad_func in zip(active_indices, active_grad_funcs):
    grad_val = grad_func(x_star)
    grad_vectors.append(grad_val)
    print(f"∇g_{i+1}(x*) = {grad_val.tolist()}")

# Step 4 & 5: Characterize and represent the normal cone
print("\n--- Step 3: Characterize the Normal Cone ---")
print("The constraint functions g_1, g_2 are convex (their Hessians are positive semidefinite),")
print("and g_3, g_4 are linear (and thus convex). Since all constraints are convex, the feasible set F is convex.")
print("\nFor a convex problem, the normal cone T_F^°(x*) is the conic hull of the gradients of the active constraints:")
print("T_F^°(x*) = { s | s = Σ_{i ∈ I(x*)} μ_i ∇g_i(x*), with μ_i ≥ 0 }")

print("\nSubstituting the calculated gradients, any vector s = (s₁, s₂, s₃) in the normal cone can be written as:")
grad1_str = f"[{grad_vectors[0][0]}, {grad_vectors[0][1]}, {grad_vectors[0][2]}]"
grad2_str = f"[{grad_vectors[1][0]}, {grad_vectors[1][1]}, {grad_vectors[1][2]}]"
grad3_str = f"[{grad_vectors[2][0]}, {grad_vectors[2][1]}, {grad_vectors[2][2]}]"
print(f"s = μ₁ * {grad1_str} + μ₂ * {grad2_str} + μ₃ * {grad3_str}, for μ₁, μ₂, μ₃ ≥ 0")

print("\n--- Step 4: Derive the Explicit Representation ---")
print("Let's analyze the components of s:")
print(f"s₁ = μ₁ * ({grad_vectors[0][0]}) + μ₂ * ({grad_vectors[1][0]}) = 2μ₁ - 2μ₂")
print(f"s₂ = μ₁ * ({grad_vectors[0][1]}) + μ₂ * ({grad_vectors[1][1]}) + μ₃ * ({grad_vectors[2][1]}) = 0")
print(f"s₃ = μ₁ * ({grad_vectors[0][2]}) + μ₂ * ({grad_vectors[1][2]}) + μ₃ * ({grad_vectors[2][2]}) = μ₃")

print("\nFrom this we deduce the properties of s₁, s₂, s₃:")
print(" - For s₁, since μ₁ ≥ 0 and μ₂ ≥ 0, the expression 2(μ₁ - μ₂) can take any real value. So, s₁ ∈ ℝ.")
print(" - For s₂, the value is always 0. So, s₂ = 0.")
print(" - For s₃, since μ₃ ≥ 0, we have s₃ ≥ 0.")

print("\nCombining these conditions, we get the explicit representation of the normal cone:")
print("T_F^°(x^*) = { s = (s₁, s₂, s₃) ∈ ℝ³ | s₁ ∈ ℝ, s₂ = 0, s₃ ≥ 0 }")
