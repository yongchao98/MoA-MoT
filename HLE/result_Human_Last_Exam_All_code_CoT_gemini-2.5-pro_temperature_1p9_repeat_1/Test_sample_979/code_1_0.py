def print_magnetic_field_solution():
    """
    This function prints the derived expressions for the magnetic field H
    inside and outside the spherical shell, corresponding to the correct answer.
    """
    
    # Define the components of the magnetic field H as strings
    H_inside_str = "H(r, θ) = (2 * μ₀ / μ) * (K₀ / (1 + (2 * μ₀ / μ))) * ẑ"
    H_outside_str = "H(r, θ) = (K₀ / (1 + (2 * μ₀ / μ))) * (R³ / r³) * (2*cos(θ) r̂ + sin(θ) θ̂)"

    # Print the solution in a formatted way
    print("The magnetic field H(r, θ) is determined in two regions:")
    print("-" * 50)
    
    print("For the region inside the sphere (0 < r < R):")
    print(H_inside_str)
    
    print("\nFor the region outside the sphere (R < r < ∞):")
    print(H_outside_str)
    print("-" * 50)
    
    # As requested, output the terms in the final equations (symbolically)
    print("\nSymbolic components of the final answer:")
    # Inside field terms
    print("H_in is proportional to a constant vector z_hat.")
    print("Coefficient for H_in: (2 * μ₀ / μ) * (K₀ / (1 + 2*μ₀/μ))")
    # Outside field terms
    print("H_out has a radial dependence of 1/r³.")
    print("Coefficient for H_out: K₀ / (1 + 2*μ₀/μ)")
    print("Angular part for H_out: (2*cos(θ)*r̂ + sin(θ)*θ̂)")


# Execute the function to display the answer
print_magnetic_field_solution()
