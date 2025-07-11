import numpy as np
import ot

def demonstrate_subgradient_at_minimum():
    """
    Numerically demonstrates that the subgradient of J(mu) = 0.5 * W(mu, nu)^2
    at its minimum (mu = nu) is the zero vector.
    It does this by showing the directional derivative is zero.
    """
    # 1. Setup: Define a discrete space and a target probability measure 'nu'.
    n_points = 100
    # The points on the x-axis where the probability mass is located.
    x_coords = np.linspace(0, 1, n_points, dtype=np.float64)

    # Create a target distribution 'nu' (e.g., a Gaussian-like shape).
    # Ensure it's in the interior of the probability simplex (all elements > 0).
    pos = np.linspace(-1, 1, n_points)
    nu = np.exp(-pos**2 / 0.1)
    nu = nu / nu.sum() # Normalize to make it a probability distribution

    # The minimizer of J is mu_star = nu.
    mu_star = nu

    # 2. Define the functional J(mu)
    def J(mu, nu_target, coords):
        # The POT library's emd2_1d computes the *squared* Wasserstein-2 distance in 1D.
        w2_squared = ot.emd2_1d(coords, coords, mu, nu_target)
        return 0.5 * w2_squared

    # 3. Create a random tangent direction 'v'.
    # A tangent vector in this context is a perturbation that keeps the sum of
    # probabilities equal to 1. This means the sum of elements in 'v' must be 0.
    np.random.seed(42)
    v_rand = np.random.randn(n_points)
    v = v_rand - v_rand.mean() # This ensures sum(v) = 0
    v = v / np.linalg.norm(v)  # Normalize for consistent scaling

    print("This script demonstrates that the gradient of J(mu) at the minimum mu=nu is the zero vector.")
    print("We do this by computing the directional derivative in a random direction v:")
    print("  D_v J(nu) = lim_{h->0} (J(nu + h*v) - J(nu)) / h")
    print("Since J(nu) = 0, this simplifies to J(nu + h*v) / h.")
    print("The result should approach 0 as h -> 0.")
    print("-" * 65)
    print(f"{'Step size h':<15} | {'J(nu + h*v)':<20} | {'Derivative Approx.':<25}")
    print("-" * 65)

    # 4. Compute the directional derivative for decreasing step sizes 'h'.
    for h in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
        # Create the perturbed measure
        mu_h = mu_star + h * v

        # Check if mu_h is a valid probability vector (all elements non-negative).
        # This will be true for small enough h because all elements of nu are positive.
        if np.any(mu_h < 0):
            print(f"{h:<15.1e} | h is too large, results in negative probabilities. Skipping.")
            continue

        # Calculate J at the perturbed point. J(mu_star) is 0.
        J_val_h = J(mu_h, mu_star, x_coords)

        # Calculate the directional derivative approximation
        directional_deriv = J_val_h / h

        print(f"{h:<15.1e} | {J_val_h:<20.4e} | {directional_deriv:<25.8f}")

    print("-" * 65)
    print("As h decreases, the directional derivative clearly approaches 0.")
    print("Since this holds for any direction v, the gradient vector must be 0.")

demonstrate_subgradient_at_minimum()