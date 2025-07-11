import numpy as np

def demonstrate_wasserstein_gradient_at_minimum():
    """
    This function demonstrates that the Wasserstein gradient of the functional
    J(mu) = 0.5 * W(mu, nu)^2 is the trivial tangent vector at its minimum.

    The demonstration uses 1D Gaussian distributions, where the gradient has a known
    analytical form.
    """

    # Define the parameters for the target distribution nu = N(m_nu, sigma_nu^2)
    m_nu = 5.0
    sigma_nu = 2.0

    print("The user's statement is: The Wasserstein regular subgradient of J(mu) = 0.5 * W(mu, nu)^2 is the trivial tangent vector at the minimum of J.")
    print("\nThis statement is TRUE. Let's demonstrate this with an example using 1D Gaussian distributions.")
    print("-" * 80)

    print(f"Let the target distribution be nu = N(mean={m_nu}, variance={sigma_nu**2}).")
    print("The functional J(mu) is minimized at mu = nu.")
    print("\nWe need to find the gradient of J at this minimum point.")

    print("\nFor a general Gaussian distribution mu = N(m, s^2), the gradient of J at mu is the vector field grad(x) given by:")
    print("grad(x) = x - T(x), where T is the optimal transport map from mu to nu.")
    print("For 1D Gaussians, T(x) = m_nu + (sigma_nu / s) * (x - m).")
    print("\nSo, the gradient formula is:")
    print("grad(x) = x - [m_nu + (sigma_nu / s) * (x - m)]")
    print("grad(x) = (m - m_nu) + (1 - sigma_nu / s) * (x - m)")
    print("-" * 80)

    print("Now, let's evaluate this gradient at the minimum, where mu = nu.")
    print(f"This means we set m = m_nu = {m_nu} and s = sigma_nu = {sigma_nu}.")

    # At the minimum, the parameters of mu are the same as nu
    m_min = m_nu
    s_min = sigma_nu

    # Calculate the coefficients of the gradient equation at the minimum
    # grad(x) = C1 + C2 * (x - m_min)
    C1 = m_min - m_nu
    C2 = 1 - sigma_nu / s_min

    print("\nSubstituting these values into the gradient formula:")
    print(f"grad(x) = ({m_min} - {m_nu}) + (1 - {sigma_nu} / {s_min}) * (x - {m_min})")
    print(f"grad(x) = ({C1}) + (1 - {sigma_nu / s_min}) * (x - {m_min})")
    print(f"grad(x) = {C1} + ({C2}) * (x - {m_min})")
    print(f"grad(x) = 0.0 + 0.0 * (x - {m_min})")
    print("grad(x) = 0")
    print("-" * 80)

    print("\nThe resulting gradient is the zero vector field, which is the trivial tangent vector in the Wasserstein space.")
    print("Since the subgradient at the minimum consists of just this gradient, the statement is confirmed to be true.")

if __name__ == '__main__':
    demonstrate_wasserstein_gradient_at_minimum()
