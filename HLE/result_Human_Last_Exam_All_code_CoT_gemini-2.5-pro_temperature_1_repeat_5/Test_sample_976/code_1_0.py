import sympy

def solve_electromagnetism_problem():
    """
    This function presents the solution for the electric potential and field
    outside a conductive sphere in a uniform electric field in a steady state.
    """
    # Define the symbolic variables used in the equations
    r, R, theta, E0, sigma1, sigma2 = sympy.symbols('r R theta E_0 sigma_1 sigma_2', real=True, positive=True)

    print("This problem asks for the electric potential Φ and electric field E outside a conductive sphere")
    print("in a uniform electric field, once a steady state has been reached.")
    print("\nThe problem is solved by using Laplace's equation (∇²Φ = 0) with the appropriate boundary conditions:")
    print("1. Continuity of potential (Φ) at the boundary r=R.")
    print("2. Continuity of the normal component of current density (J = σE) at r=R.")
    print("3. The potential far from the sphere approaches that of the uniform field (-E₀*r*cos(θ)).")
    print("\nThe derivation yields the expressions in choice B. We will now display the results for the region outside the sphere (r > R).")

    # --- Electric Potential Φ for r > R ---
    print("\n-------------------------------------------------------------")
    print("Electric Potential Φ(r, θ) for r > R")
    print("-------------------------------------------------------------")
    
    # The potential is composed of the external field potential and a dipole term.
    # The coefficient of the dipole term is derived from the boundary conditions.
    # Φ(r, θ) = -E₀ * ( r - ( (σ₁ - σ₂) / (σ₁ + 2*σ₂) ) * R³ / r² ) * cos(θ)
    
    print("The final equation for the potential outside the sphere is:")
    print("Φ(r, θ) = -E_0 * ( r - term_coefficient * (R**3 / r**2) ) * cos(θ)")
    print("where the term coefficient is (σ₁ - σ₂) / (σ₁ + 2*σ₂).")
    print("\nFull equation:")
    print("Φ(r, θ) = -E_0 * (r - ((sigma_1 - sigma_2) * R**3) / ((sigma_1 + 2*sigma_2) * r**2)) * cos(theta)")


    # --- Electric Field E for r > R ---
    print("\n-------------------------------------------------------------")
    print("Electric Field E(r, θ) for r > R")
    print("-------------------------------------------------------------")
    print("The electric field is found by taking the negative gradient of the potential, E = -∇Φ.")
    print("The field has a radial component (E_r) and an angular component (E_θ).")
    
    # Radial Component E_r
    print("\nRadial Component (E_r):")
    print("E_r(r, θ) = E_0 * (1 + 2 * term_coefficient * (R**3 / r**3)) * cos(θ)")
    print("where the term coefficient is (σ₁ - σ₂) / (σ₁ + 2*σ₂).")
    print("\nFull equation:")
    print("E_r(r, θ) = E_0 * (1 + (2 * (sigma_1 - sigma_2) * R**3) / ((sigma_1 + 2*sigma_2) * r**3)) * cos(theta)")

    # Angular Component E_θ
    print("\nAngular Component (E_θ):")
    print("E_θ(r, θ) = -E_0 * (1 - term_coefficient * (R**3 / r**3)) * sin(θ)")
    print("where the term coefficient is (σ₁ - σ₂) / (σ₁ + 2*σ₂).")
    print("\nFull equation:")
    print("E_θ(r, θ) = -E_0 * (1 - ((sigma_1 - sigma_2) * R**3) / ((sigma_1 + 2*sigma_2) * r**3)) * sin(theta)")

    print("\n\nThese results for the potential and field components outside the sphere match answer choice B.")

if __name__ == '__main__':
    solve_electromagnetism_problem()