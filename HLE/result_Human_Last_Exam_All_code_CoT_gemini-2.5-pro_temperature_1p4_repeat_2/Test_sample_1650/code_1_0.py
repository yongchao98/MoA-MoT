import sympy
from sympy import exp, integrate, oo, pi, Symbol, simplify

def solve_overlap_integral():
    """
    Calculates the overlap integral for two 2s hydrogenic orbitals in H2+
    using symbolic integration with SymPy.
    """
    # Define symbolic variables
    zeta = Symbol('ζ', positive=True, real=True)  # Effective nuclear charge
    R = Symbol('R', positive=True, real=True)      # Internuclear distance
    lmbda = Symbol('λ', real=True)
    mu = Symbol('μ', real=True)

    # Use a single composite variable rho = ζ*R for the final expression
    rho = Symbol('ρ', positive=True, real=True)

    # The 2s hydrogenic wavefunction is ψ_2s = N * (2 - ζr) * exp(-ζr/2)
    # The squared normalization constant is N^2 = ζ^3 / (32π)
    N_squared = zeta**3 / (32 * pi)

    # In elliptical coordinates, the product of the polynomial parts is:
    # [2 - ζ*r_a] * [2 - ζ*r_b] = 4 - 2*ζ*R*λ + (ζ^2*R^2/4)*(λ^2 - μ^2)
    # We substitute ζR with ρ.
    poly_part = 4 - 2*rho*lmbda + (rho**2 / 4) * (lmbda**2 - mu**2)

    # The product of the exponential parts is:
    # exp(-ζ*r_a/2) * exp(-ζ*r_b/2) = exp(-ζ*R*λ/2)
    exp_part = exp(-rho * lmbda / 2)

    # The volume element dτ in elliptical coordinates is (R^3/8) * (λ^2 - μ^2) dλ dμ dφ.
    # We substitute R = ρ/ζ.
    d_tau_factor = ((rho/zeta)**3 / 8) * (lmbda**2 - mu**2)

    # The full integrand for the overlap integral S = ∫ N^2 * (ψ_a*ψ_b) dτ
    integrand = N_squared * poly_part * exp_part * d_tau_factor

    # The integral is over φ from 0 to 2π, μ from -1 to 1, and λ from 1 to ∞.

    # Step 1: Integrate with respect to φ from 0 to 2π
    # Since the integrand is independent of φ, this is a multiplication by 2π.
    integral_after_phi = (2 * pi * integrand).simplify()

    # Step 2: Integrate with respect to μ from -1 to 1
    integral_after_mu = integrate(integral_after_phi, (mu, -1, 1))

    # Step 3: Integrate with respect to λ from 1 to infinity
    S = integrate(integral_after_mu, (lmbda, 1, oo))
    
    # Simplify the final expression
    S_simplified = simplify(S)

    # Print the final result in a clear, readable format.
    # The output shows the final equation with all its terms and numbers.
    final_expression_str = str(S_simplified).replace('rho', 'ρ')
    print("The analytical expression for the overlap integral S(2s, 2s) is:")
    print(f"S(ρ) = {final_expression_str}")
    print("\nwhere ρ = ζR.")

if __name__ == '__main__':
    solve_overlap_integral()
