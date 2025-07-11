import sympy as sp

def solve_electrosatics_problem():
    """
    This function defines and prints the expressions for the electric potential and
    electric field outside a conductive sphere in a uniform electric field,
    in a steady state, using symbolic mathematics.
    """
    # Define symbolic variables
    E0, r, R, theta, sigma1, sigma2 = sp.symbols('E_0 r R theta sigma_1 sigma_2')
    
    # --- Potential and Field outside the sphere (r > R) ---

    # Electric Potential (Phi)
    # Φ(r, θ) = -E₀ * (r - (σ₁ - σ₂) * R³ / ((σ₁ + 2*σ₂) * r²)) * cos(θ)
    phi_out_factor = E0 * (r - (sigma1 - sigma2) * R**3 / ((sigma1 + 2 * sigma2) * r**2))
    phi_out = -phi_out_factor * sp.cos(theta)

    # Electric Field (E = -∇Φ)
    # E_r = -∂Φ/∂r
    E_out_r = -sp.diff(phi_out, r)
    
    # E_θ = -(1/r) * ∂Φ/∂θ
    E_out_theta = -(1/r) * sp.diff(phi_out, theta)

    # Print the results for the region outside the sphere (r > R)
    print("The electric potential and electric field in the region outside the sphere (r > R) are:")
    
    print("\nElectric Potential Φ(r, θ):")
    sp.pprint(phi_out, use_unicode=True)
    
    print("\nElectric Field E(r, θ) = E_r r_hat + E_θ θ_hat:")
    print("\nRadial component (E_r):")
    # The default simplified form from sympy is a bit different from the answer choices,
    # but mathematically equivalent. We'll build the expression as in the choices for clarity.
    E_r_expr = E0 * (1 + 2 * (sigma1 - sigma2) * R**3 / ((sigma1 + 2*sigma2) * r**3)) * sp.cos(theta)
    sp.pprint(E_r_expr, use_unicode=True)

    print("\nAngular component (E_θ):")
    E_theta_expr = -E0 * (1 - (sigma1 - sigma2) * R**3 / ((sigma1 + 2*sigma2) * r**3)) * sp.sin(theta)
    sp.pprint(E_theta_expr, use_unicode=True)

    print("\nThese results match the expressions for r > R in answer choice B.")

solve_electrosatics_problem()