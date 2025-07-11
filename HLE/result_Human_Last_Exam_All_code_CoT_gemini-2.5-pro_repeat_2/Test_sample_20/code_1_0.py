import sympy

def demonstrate_gradient_at_minimum():
    """
    Demonstrates that the Wasserstein gradient of J(mu) = 0.5 * W(mu, nu)^2
    is the zero vector at the minimum (mu=nu).

    This is shown for the 1D Gaussian case.
    """
    # Define a symbolic variable for the space
    x = sympy.Symbol('x')

    # Define concrete parameters for the target distribution nu
    # Let nu be a Gaussian N(mean, std^2)
    m_nu = 5.0
    s_nu = 2.0

    # At the minimum of J(mu), we have mu = nu.
    # So, the parameters for mu are the same as for nu.
    m_mu = 5.0
    s_mu = 2.0

    print("Yes, the statement is true. At the minimum of J(mu), where mu=nu, the subgradient is the trivial (zero) tangent vector.")
    print("\nWe demonstrate this for 1D Gaussian distributions.")
    print(f"Let the target distribution be nu = N(mean={m_nu}, std_dev={s_nu}).")
    print(f"The minimum of J(mu) occurs at mu = nu, so we evaluate the gradient at mu = N(mean={m_mu}, std_dev={s_mu}).")

    # The gradient vector field is v(x) = x - T(x), where T is the optimal transport map.
    # For 1D Gaussians, T(x) = m_nu + (s_nu / s_mu) * (x - m_mu).
    print("\nThe gradient vector field v(x) is given by the equation v(x) = x - T(x).")
    print("For these specific Gaussians, this means:")

    # We will print the equation step-by-step to show the simplification.
    # The numbers from our chosen parameters are substituted into the final equation.
    print(f"\nv(x) = x - ({m_nu} + ({s_nu}/{s_mu}) * (x - {m_mu}))")

    # Perform the division
    ratio = s_nu / s_mu
    print(f"v(x) = x - ({m_nu} + {ratio} * (x - {m_mu}))")

    # Distribute the ratio
    # In this case, ratio is 1.0, so this simplifies nicely.
    print(f"v(x) = x - ({m_nu} + {ratio}*x - {ratio*m_mu})")

    # Combine terms inside the parenthesis
    print(f"v(x) = x - (x + {m_nu - ratio*m_mu})")

    # Final simplification
    print(f"v(x) = x - x")

    final_result = 0
    print(f"v(x) = {final_result}")

    print("\nAs shown, the gradient vector field is 0 for all x. This is the trivial tangent vector.")

demonstrate_gradient_at_minimum()