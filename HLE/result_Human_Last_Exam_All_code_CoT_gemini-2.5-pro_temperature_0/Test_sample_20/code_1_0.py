import numpy as np
import ot

# This script demonstrates that the Wasserstein gradient of the functional
# J(mu) = 0.5 * W(mu, nu)^2 is the zero tangent vector at its minimum (mu = nu).

# 1. Define a discrete probability measure nu.
# We create a measure supported on N points in d-dimensional space.
N_POINTS = 10
DIM = 2
np.random.seed(42) # for reproducibility

# The support points of the measure nu
support_nu = np.random.rand(N_POINTS, DIM)
# The weights of the measure nu (uniform)
weights_nu = np.ones(N_POINTS) / N_POINTS

# 2. Set mu = nu to evaluate the gradient at the minimum of J.
support_mu = support_nu
weights_mu = weights_nu

print("--- Setup ---")
print(f"The measure nu is defined on {N_POINTS} points in {DIM}D space.")
print("We are calculating the Wasserstein gradient of J(mu) at its minimum, where mu = nu.")
print("\nSupport points (x_i) of the measure mu=nu:")
print(support_mu)

# 3. Compute the optimal transport map T from mu to nu.
# The cost matrix C_ij = ||x_i - y_j||^2. Since mu=nu, x=y.
cost_matrix = ot.dist(support_mu, support_nu, metric='sqeuclidean')

# Solve the Earth Mover's Distance problem to get the optimal transport plan (gamma).
# For mu=nu, the cost to transport mass from x_i to x_i is 0, so the optimal
# plan should be a diagonal matrix where gamma_ii = weights_mu[i].
ot_plan = ot.emd(weights_mu, weights_nu, cost_matrix)

# The optimal transport map T sends each point x_i in mu's support to T(x_i).
# T(x_i) is calculated as a weighted average of points in nu's support:
# T(x_i) = sum_j (ot_plan[i,j] / weights_mu[i]) * support_nu[j]
# Since ot_plan is diagonal, T(x_i) should be equal to x_i.
transported_points = (ot_plan / weights_mu[:, np.newaxis]) @ support_nu

print("\n--- Calculation ---")
print("The gradient vector field is v(x) = x - T(x).")
print("We compute this for each support point x_i:")

# 4. Compute the gradient vector field v_i = x_i - T(x_i) and print the equation.
gradient_field = support_mu - transported_points

is_trivial = True
for i in range(N_POINTS):
    x_i = support_mu[i]
    Tx_i = transported_points[i]
    v_i = gradient_field[i]
    
    # Check if the resulting vector is close to zero
    if not np.allclose(v_i, 0):
        is_trivial = False
        
    print(f"\nPoint {i}: v_{i} = x_{i} - T(x_{i})")
    # Format numbers for cleaner printing
    x_i_str = np.array2string(x_i, formatter={'float_kind':lambda x: "%.4f" % x})
    Tx_i_str = np.array2string(Tx_i, formatter={'float_kind':lambda x: "%.4f" % x})
    v_i_str = np.array2string(v_i, formatter={'float_kind':lambda x: "%.4f" % x})
    print(f"         v_{i} = {x_i_str} - {Tx_i_str} = {v_i_str}")

# 5. Final conclusion based on the numerical result.
print("\n--- Conclusion ---")
if is_trivial:
    print("The gradient vector field is numerically zero at all support points.")
    print("This confirms that the Wasserstein regular subgradient at the minimum is the trivial (zero) tangent vector.")
    print("\nTherefore, the statement is TRUE.")
else:
    print("The gradient vector field is not zero.")
    print("Therefore, the statement is FALSE.")
