import numpy as np
import ot

def demonstrate_wasserstein_gradient_at_minimum():
    """
    This function numerically demonstrates that the Wasserstein gradient of
    J(mu) = 0.5 * W(mu, nu)^2 is the zero vector at the minimum (mu = nu).
    """
    # Number of points in the discrete distributions
    n = 5
    
    # Define the target probability measure nu
    # Support points of nu
    y = np.linspace(0, 10, n).reshape((n, 1))
    # Weights of nu (uniform)
    b = np.ones((n,)) / n

    # The minimum of J is at mu = nu.
    # So we define mu to be the same as nu.
    x = y
    a = b

    print("This script demonstrates that the Wasserstein gradient of the squared")
    print("Wasserstein distance J(mu) = 0.5 * W(mu, nu)^2 is the trivial (zero) vector")
    print("at its minimum, which occurs when mu = nu.")
    print("-" * 50)
    
    # Compute the cost matrix M_{i,j} = ||x_i - y_j||^2
    M = ot.dist(x, y, metric='sqeuclidean')

    # Compute the optimal transport coupling matrix (plan) gamma
    # For discrete measures, this is equivalent to solving the Earth Mover's Distance problem.
    gamma = ot.emd(a, b, M)

    # From the coupling matrix, we can find the optimal transport map T.
    # For a source point x_i, the mapped point T(x_i) is the weighted average of target points,
    # T(x_i) = (1/a_i) * sum_j(gamma_ij * y_j)
    # Since a_i = 1/n and gamma is a permutation matrix scaled by 1/n, this simplifies.
    # T_x will be a matrix where T_x[i] = T(x_i).
    T_x = n * gamma.dot(y)

    # The Wasserstein gradient vector field at the points of mu is given by x - T(x)
    gradient_vector = x - T_x

    # --- Output the results ---
    # The "final equation" is the computed zero vector.
    print(f"Support points of the target measure nu, y:\n{y.flatten()}\n")
    print(f"Support points of the measure mu at the minimum (mu=nu), x:\n{x.flatten()}\n")
    print("The optimal transport map T sends each point in mu's support to the corresponding point in nu's support.")
    print(f"Mapped points T(x):\n{T_x.flatten()}\n")
    print("The gradient vector at the minimum is calculated as x - T(x):")
    print(f"Gradient vector = {x.flatten()} - {T_x.flatten()}")
    print("Final Result:")
    # Print each number in the final equation
    print("[", end="")
    for i, val in enumerate(gradient_vector.flatten()):
        print(f"{val:.1f}", end="")
        if i < len(gradient_vector.flatten()) - 1:
            print(", ", end="")
    print("]\n")
    print("As shown, the gradient at the minimum is the trivial (zero) tangent vector.")

if __name__ == '__main__':
    demonstrate_wasserstein_gradient_at_minimum()
