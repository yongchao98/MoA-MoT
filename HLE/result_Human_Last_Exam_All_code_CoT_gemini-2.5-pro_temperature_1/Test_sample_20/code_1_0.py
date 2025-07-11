import numpy as np

def illustrate_wasserstein_gradient():
    """
    Illustrates that the Wasserstein gradient of J(mu) = 0.5 * W_2(mu, nu)^2
    approaches the trivial (zero) tangent vector as mu approaches nu.

    We use a simple 1D Gaussian case where both mu and nu are normal distributions.
    The gradient has a simple analytical form in this case.
    """

    # Let nu be the target distribution, a standard normal N(0, 1).
    mean_nu = 0.0
    var_nu = 1.0

    print("This script demonstrates the behavior of the Wasserstein gradient for the functional:")
    print("J(mu) = (1/2) * W(mu, nu)^2")
    print(f"The minimum of J is at mu = nu, where nu is the Gaussian N({mean_nu}, {var_nu}).")
    print("\nWe will consider a family of distributions mu_a = N(a, 1) which approach nu as 'a' -> 0.")
    print("-" * 70)
    print("Theoretical Background for 1D Gaussians:")
    print("  mu = N(m1, s^2), nu = N(m2, s^2)  (note: same variance for simplicity)")
    print("  Squared Wasserstein distance: W(mu, nu)^2 = (m1 - m2)^2")
    print("  Optimal Transport Map from mu to nu: T(x) = x - (m1 - m2)")
    print("  Wasserstein Gradient Field of J(mu): grad(J)(x) = x - T(x) = m1 - m2")
    print("-" * 70)

    print("\nLet's observe the gradient as the mean 'a' of mu_a approaches 0:")
    print("\n{:<15} {:<25} {:<30}".format("Mean 'a' of mu_a", "Value of J(mu_a)", "Gradient Field (a constant value)"))

    # Consider a sequence of 'a' values (means of mu_a) approaching 0.
    a_values = [1.0, 0.5, 0.1, 0.01, 0.001, 0.0]

    for a in a_values:
        # mu_a is the distribution N(a, 1).
        mean_mu = a
        # For simplicity, we keep the variance the same. The logic holds for varying variance too.
        var_mu = 1.0

        # The squared Wasserstein distance is W^2 = (m1 - m2)^2 + (s1 - s2)^2
        # Here s1=s2, so W^2 = (m1-m2)^2
        w2_squared = (mean_mu - mean_nu)**2

        # The value of the functional J(mu) = 0.5 * W^2
        j_value = 0.5 * w2_squared

        # The gradient is the constant vector field with value (mean_mu - mean_nu)
        gradient_value = mean_mu - mean_nu

        print("{:<15.4f} {:<25.8f} {:<30.4f}".format(a, j_value, gradient_value))

    print("\nAs 'a' -> 0, mu_a -> nu. We observe that J(mu_a) -> 0 and the gradient -> 0.")
    print("At the minimum (a=0), the gradient is exactly 0, the trivial tangent vector.")

if __name__ == '__main__':
    illustrate_wasserstein_gradient()