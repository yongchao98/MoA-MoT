def print_solution():
    """
    This function prints the derived expressions for the electric potential and
    electric field outside the sphere (r > R).
    """

    # Define the expressions as strings for clear printing.
    # Note: σ₁ is sigma1, σ₂ is sigma2, E₀ is E0, R is the radius.

    # Potential outside the sphere (r > R)
    phi_out_str = "Φ(r, θ) = -E0 * (r - (σ1 - σ2) * R^3 / ((σ1 + 2*σ2) * r^2)) * cos(θ)"

    # Electric field outside the sphere (r > R)
    # E_r component
    E_r_out_str = "E_r(r, θ) = E0 * [1 + 2*(σ1 - σ2) * R^3 / ((σ1 + 2*σ2) * r^3)] * cos(θ)"
    # E_theta component
    E_theta_out_str = "E_θ(r, θ) = -E0 * [1 - (σ1 - σ2) * R^3 / ((σ1 + 2*σ2) * r^3)] * sin(θ)"
    E_vec_out_str = "E_vector(r, θ) = E_r(r, θ) * r_hat + E_θ(r, θ) * θ_hat"

    print("The solution for the electric potential and field outside the sphere (r > R) is as follows:")
    
    print("\n--- Electric Potential Φ(r, θ) ---")
    print(phi_out_str)

    print("\n--- Electric Field E(r, θ) ---")
    print("The electric field is a vector with radial (r_hat) and polar (θ_hat) components:")
    print(E_vec_out_str)
    
    print("\nRadial component:")
    print(E_r_out_str)
    
    print("\nPolar component:")
    print(E_theta_out_str)
    
    print("\nThese expressions correspond to the solution given in Answer Choice B.")

# Execute the function to print the solution
print_solution()