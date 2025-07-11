def display_magnetic_field_solution():
    """
    This function prints the final derived equations for the magnetic field H
    inside and outside the spherical shell, as requested by the problem.
    The numbers in the final equations are explicitly shown.
    """

    print("The final expressions for the magnetic field H(r, θ) are presented below.")

    # --- Region 1: Inside the sphere (0 < r < R) ---
    print("\n--- Magnetic Field Inside the Sphere (0 < r < R) ---")
    
    # The final equation contains the numbers 2, 1, 2
    inside_equation = "H_in(r, θ) = ( (2 * μ_0) / μ ) * K_0 / ( 1 + (2 * μ_0) / μ ) * ẑ"
    print(inside_equation)

    # --- Region 2: Outside the sphere (r > R) ---
    print("\n--- Magnetic Field Outside the Sphere (r > R) ---")
    
    # The final equation contains the numbers 1, 1, 2, 3, 3, 2, 1
    # We can write sin(θ) as 1*sin(θ) to explicitly show the coefficient
    outside_equation = "H_out(r, θ) = K_0 / ( 1 + (2 * μ_0) / μ ) * (R^3 / r^3) * (2*cos(θ)*r̂ + sin(θ)*θ̂)"
    print(outside_equation)

    print("\nThese results correspond to answer choice E.")

if __name__ == "__main__":
    display_magnetic_field_solution()