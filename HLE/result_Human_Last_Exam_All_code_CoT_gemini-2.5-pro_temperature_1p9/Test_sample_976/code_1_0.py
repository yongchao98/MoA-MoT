def print_solution():
    """
    This function prints the final expressions for the electric potential and
    electric field outside a conductive sphere in a uniform electric field
    in a steady state.
    """

    # The problem is a classic electrostatics boundary value problem.
    # The governing equation is Laplace's equation,nabla^2(Phi) = 0, for the
    # electric potential Phi in each region (inside and outside the sphere).
    # The key boundary condition at the sphere's surface (r=R) for this
    # steady-state current problem is the continuity of the normal component of
    # the current density (J_n = sigma * E_n), which means:
    # sigma_1 * d(Phi_in)/dr = sigma_2 * d(Phi_out)/dr at r=R.
    # The solution is analogous to that of a dielectric sphere in a uniform
    # field, with permittivities (epsilon_1, epsilon_2) replaced by
    # conductivities (sigma_1, sigma_2).
    #
    # The correct derivation leads to the formulas presented in option B.
    # This script will display the formulas for the region outside the sphere (r > R).

    print("The electric potential and electric field in the region outside the sphere (r > R) are:")
    print("-" * 80)

    # --- Electric Potential Outside (r > R) ---
    phi_out = "Phi(r, theta) = -E_0 * (r - (sigma_1 - sigma_2)*R**3 / ((sigma_1 + 2*sigma_2)*r**2)) * cos(theta)"
    print("Electric Potential:")
    print(phi_out)
    print("") # for spacing

    # --- Electric Field Outside (r > R) ---
    # The electric field is E = -grad(Phi). It has a radial (r) and a polar (theta) component.
    e_r_out = "E_r = E_0 * (1 + 2*(sigma_1 - sigma_2)*R**3 / ((sigma_1 + 2*sigma_2)*r**3)) * cos(theta)"
    e_theta_out = "E_theta = -E_0 * (1 - (sigma_1 - sigma_2)*R**3 / ((sigma_1 + 2*sigma_2)*r**3)) * sin(theta)"

    print("Electric Field (vector components):")
    print(e_r_out)
    print(e_theta_out)
    print("\nWhere E_vector = E_r * r_hat + E_theta * theta_hat")
    print("-" * 80)

if __name__ == "__main__":
    print_solution()