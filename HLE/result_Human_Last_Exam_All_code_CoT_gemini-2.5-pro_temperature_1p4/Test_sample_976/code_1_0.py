def print_solution():
    """
    This function prints the final derived expressions for the potential and electric field
    outside the sphere, as requested by the user.
    """
    
    # Final Potential outside the sphere (r > R)
    phi_out = r"Φ(r, θ) = -E₀ ( r - ( (σ₁ - σ₂) R³ ) / ( (σ₁ + 2σ₂)r² ) ) cos(θ)"
    
    # Final Electric Field outside the sphere (r > R)
    E_r_out = r"E_r = E₀ [ 1 + ( 2(σ₁ - σ₂) R³ ) / ( (σ₁ + 2σ₂) r³ ) ] cos(θ)"
    E_theta_out = r"E_θ = -E₀ [ 1 - ( (σ₁ - σ₂) R³ ) / ( (σ₁ + 2σ₂) r³ ) ] sin(θ)"
    E_vec_out = r"vec(E)(r, θ) = E_r r_hat + E_θ θ_hat"

    print("The final expressions for the potential and electric field in the region outside the sphere (r > R) are:")
    print("-" * 80)
    print("Potential Φ(r,θ) for r > R:")
    print(phi_out)
    print("\nElectric Field E(r,θ) for r > R:")
    print(E_vec_out)
    print("where:")
    print(f"  {E_r_out}")
    print(f"  {E_theta_out}")
    print("-" * 80)
    print("These expressions match the results given in answer choice B.")

print_solution()