import numpy as np
from scipy.integrate import quad

def solve_critical_chemical_potential():
    """
    Calculates the critical chemical potential for a holographic superconductor
    in 5D Einstein-Gauss-Bonnet gravity using a variational method.
    """
    
    # Model Parameters
    lambda_gb = 0.1
    m_sq = -3  # Corresponds to a dimension 3 operator like \bar{q}q

    # The EGB metric function f(u) in dimensionless coordinates u=z/z_h
    # where z_h is the horizon position.
    def f(u):
        # Using a small epsilon to avoid division by zero at u=0 for numerical stability,
        # though the integrand is well-behaved.
        eps = 1e-9
        u = np.maximum(u, eps)
        # We define a = 1 - 4*\lambda_gb. The argument of the square root is a + 4*\lambda_gb*u^4
        a = 1 - 4 * lambda_gb
        return (1 / (2 * lambda_gb)) * (1 - np.sqrt(a + 4 * lambda_gb * u**4))

    # We use a trial function F(u) = u^3(1-u) which respects the boundary conditions
    # for a Delta=3 operator (F(u) ~ u^3 at boundary u=0 and F(1)=0 at horizon).
    F = lambda u: u**3 * (1 - u)
    dF_du = lambda u: 3 * u**2 - 4 * u**3
    
    # Numerator of the Sturm-Liouville variational formula for (\mu_c z_h)^2
    def numerator_integrand(u):
        val = (f(u) / u**2) * dF_du(u)**2 + (m_sq / u**4) * F(u)**2
        return val

    # Denominator of the Sturm-Liouville variational formula
    def denominator_integrand(u):
        val = ((1 - u**2)**2 / (u**2 * f(u))) * F(u)**2
        return val

    # Perform the numerical integrations
    # We integrate from a small epsilon to 1 to avoid any potential numerical issues at u=0.
    eps = 1e-9
    numerator_val, _ = quad(numerator_integrand, eps, 1)
    denominator_val, _ = quad(denominator_integrand, eps, 1)

    # The eigenvalue from the variational principle is \Lambda = (\mu_c z_h)^2
    lambda_eigenvalue = numerator_val / denominator_val
    
    # The ratio \mu_c / T_c can be calculated from the eigenvalue.
    # The temperature is T_c = |f'(1)| / (4\pi z_h) = 1/(\pi z_h) for this model.
    # Therefore, \mu_c/T_c = \pi * sqrt(\Lambda)
    mu_c_over_Tc = np.pi * np.sqrt(lambda_eigenvalue)
    
    # Print the result and the equation
    print("In the context of bottom-up holographic models, the critical chemical potential is often computed as a dimensionless ratio with respect to the critical temperature, mu_c / T_c.")
    print("This can be estimated using the Sturm-Liouville variational method with a trial function F(u).")
    print("The formula for the ratio squared is:")
    print("(mu_c / T_c)^2 = pi^2 * [Integral_0^1( Numerator )du] / [Integral_0^1( Denominator )du]\n")
    print("Where:")
    print("Numerator = (f(u)/u^2)(dF/du)^2 + (m^2/u^4)F(u)^2")
    print("Denominator = ((1-u^2)^2 / (u^2 f(u))) F(u)^2")
    print(f"With parameters m^2 = {m_sq}, lambda_GB = {lambda_gb} and trial function F(u) = u^3(1-u).\n")

    print("--- Calculation ---")
    print(f"Value of the numerator integral: {numerator_val:.4f}")
    print(f"Value of the denominator integral: {denominator_val:.4f}")
    print(f"Eigenvalue Lambda = (mu_c * z_h)^2 = {numerator_val:.4f} / {denominator_val:.4f} = {lambda_eigenvalue:.4f}")
    print(f"The critical chemical potential over critical temperature, mu_c/T_c = pi * sqrt(Lambda) = 3.14159 * sqrt({lambda_eigenvalue:.4f}) = {mu_c_over_Tc:.4f}")
    print("\n--- Final Answer ---")
    print(f"The value of the critical chemical potential to temperature ratio is: {mu_c_over_Tc:.4f}")
    
    return mu_c_over_Tc

# Run the calculation and store the final answer
final_answer = solve_critical_chemical_potential()