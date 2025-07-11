def display_final_equation():
    """
    Constructs and prints the correct symbolic formula for the magnetic field,
    which corresponds to answer choice E.
    """
    
    # Common coefficient part in both regional solutions
    coefficient = "K_0 / (1 + (2 * mu_0 / mu))"
    
    # Field inside the sphere (r < R)
    inside_factor = "(2 * mu_0 / mu)"
    inside_field = f"    {inside_factor} * {coefficient} * z_hat"
    
    # Field outside the sphere (r > R)
    outside_spatial_part = "(R^3 / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)"
    outside_field = f"    {coefficient} * {outside_spatial_part}"
    
    # Print the full piecewise solution in a readable format
    print("The final solution for the magnetic field H(r, theta) is:")
    print("H(r, theta) = ")
    print("{")
    print(f"{inside_field}      (for 0 < r < R)")
    print("")
    print(f"{outside_field}      (for R < r < infinity)")
    print("}")
    print("\nThis mathematical expression matches answer choice E.")

# Execute the function to display the result.
display_final_equation()