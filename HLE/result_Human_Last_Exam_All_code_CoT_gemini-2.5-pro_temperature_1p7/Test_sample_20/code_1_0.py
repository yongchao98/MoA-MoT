import numpy as np
import ot

# This script illustrates that the subgradient of the squared Wasserstein distance
# functional J(mu) = 0.5 * W(mu, nu)^2 is the trivial tangent vector at its minimum.
# The minimum is achieved when mu = nu. We will show that the gradient vector field at this point is zero.

# 1. Define a reference probability measure 'nu' in R^2.
# We will use a discrete measure for simplicity.
# Let's define 4 points in 2D space.
nu_locations = np.array([
    [0., 0.],
    [1., 3.],
    [4., 2.],
    [5., 5.]
])
# Let's assign uniform weights to these points.
n_points = nu_locations.shape[0]
nu_weights = np.ones(n_points) / n_points

print("Reference measure nu:")
print("Locations:\n", nu_locations)
print("Weights:\n", nu_weights)
print("-" * 30)

# The functional J(mu) is minimized at mu = nu. We want to find the subgradient at this point.
# Let's set mu = nu.
mu_locations = nu_locations
mu_weights = nu_weights

# 2. Calculate the cost matrix.
# The cost is the squared Euclidean distance between points.
# M[i, j] = ||mu_locations[i] - nu_locations[j]||^2
M = ot.dist(mu_locations, nu_locations, metric='sqeuclidean')

# 3. Compute the optimal transport plan (coupling).
# Since mu = nu, we are transporting nu to itself.
optimal_plan = ot.emd(mu_weights, nu_weights, M)
W_sq = np.sum(optimal_plan * M)

print("Optimal transport plan from nu to nu:\n", optimal_plan)
print(f"Squared Wasserstein distance W(nu, nu)^2: {W_sq:.4f}")
print("-" * 30)

# 4. Analyze the transport plan and compute the gradient.
# The subgradient is given by vectors v(x) = x - T(x), where T is the optimal transport map.
# For discrete measures, the transport plan P_ij gives the mass transported from mu_i to nu_j.
# The optimal plan should be a diagonal matrix (scaled by weights), corresponding to the identity map T(x_i) = x_i.
is_identity = np.allclose(optimal_plan, np.diag(mu_weights))
print(f"Is the optimal plan the identity plan? {is_identity}")

# Let's compute the gradient vector field v at each point of mu's support.
# The gradient vector at mu_locations[i] is mu_locations[i] - T(mu_locations[i]).
# We can compute T(mu_locations[i]) for each i from the optimal plan.
transported_locations = np.zeros_like(mu_locations)
for i in range(n_points):
    # For each source point i, T(mu_i) is the weighted average of target points,
    # where weights are from the i-th row of the transport plan, normalized by the mass at i.
    if mu_weights[i] > 1e-9: # Check for non-zero mass to avoid division by zero
        transported_locations[i] = optimal_plan[i, :] @ nu_locations / mu_weights[i]

gradient_field = mu_locations - transported_locations

print("Transported locations T(x_i) (should be identical to x_i):\n", transported_locations)
print("Resulting gradient field v(x_i) = x_i - T(x_i):\n", gradient_field)
print("-" * 30)

# Check if the gradient field is the zero vector field (trivial tangent vector).
is_trivial = np.allclose(gradient_field, 0)
print(f"Is the gradient field the trivial tangent vector? {is_trivial}")
print("\nConclusion: The subgradient at the minimum (mu=nu) is the trivial tangent vector.")