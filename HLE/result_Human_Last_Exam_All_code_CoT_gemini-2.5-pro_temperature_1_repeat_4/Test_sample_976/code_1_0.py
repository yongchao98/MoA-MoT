def print_final_equations():
    """
    This function prints the final equations for the electric potential and electric field
    in the region outside the sphere (r > R), corresponding to the correct answer.
    """
    
    # Header for the output
    print("The final expressions for the electric potential Phi and electric field E")
    print("in the region outside the sphere (r > R) are:")
    print("-" * 60)

    # Define the equations as strings for clear formatting
    # Note: sigma_1 and sigma_2 refer to the conductivities.
    # The dielectric permittivities epsilon_1 and epsilon_2 do not appear in the
    # steady-state solution.
    
    potential_outside = "Phi(r, theta) = -E_0 * (r - ((sigma_1 - sigma_2) * R**3) / ((sigma_1 + 2*sigma_2) * r**2)) * cos(theta)"
    
    field_radial_component = "E_r(r, theta) = E_0 * (1 + (2 * (sigma_1 - sigma_2) * R**3) / ((sigma_1 + 2*sigma_2) * r**3)) * cos(theta)"
    
    field_polar_component = "E_theta(r, theta) = -E_0 * (1 - ((sigma_1 - sigma_2) * R**3) / ((sigma_1 + 2*sigma_2) * r**3)) * sin(theta)"

    # Print the Electric Potential
    print("\nElectric Potential (r > R):")
    print(potential_outside)

    # Print the Electric Field
    print("\nElectric Field (r > R) in spherical components (E = E_r * r_hat + E_theta * theta_hat):")
    print("\n[Radial component E_r]")
    print(field_radial_component)
    
    print("\n[Polar component E_theta]")
    print(field_polar_component)
    
    print("\n" + "-" * 60)
    print("These expressions correspond to Answer Choice B.")

# Execute the function to display the answer
print_final_equations()