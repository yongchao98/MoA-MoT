def print_solution():
    """
    This function prints the final expressions for the electric potential and
    electric field outside the sphere (r > R).
    """

    # Define the variables for clarity in the output strings
    E0 = "E_0"
    r = "r"
    R = "R"
    theta = "theta"
    sigma1 = "sigma_1"
    sigma2 = "sigma_2"

    # Construct the string for the potential Phi
    potential_str = (f"Phi(r, theta) = -{E0} * ( {r} - (({sigma1} - {sigma2}) * {R}**3) /"
                     f" (({sigma1} + 2*{sigma2}) * {r}**2) ) * cos({theta})")

    # Construct the strings for the components of the Electric Field E
    # Radial component
    E_r_str = (f"E_r = {E0} * [ 1 + (2*({sigma1} - {sigma2}) * {R}**3) /"
               f" (({sigma1} + 2*{sigma2}) * {r}**3) ] * cos({theta})")

    # Tangential component
    E_theta_str = (f"E_theta = -{E0} * [ 1 - (({sigma1} - {sigma2}) * {R}**3) /"
                   f" (({sigma1} + 2*{sigma2}) * {r}**3) ] * sin({theta})")

    print("The final expressions for the region outside the sphere (r > R) are:")
    print("\nElectric Potential:")
    print(potential_str)
    print("\nElectric Field (E = E_r * r_hat + E_theta * theta_hat):")
    print("Radial Component:")
    print(E_r_str)
    print("\nTangential Component:")
    print(E_theta_str)
    print("\nThese expressions correspond to the solution given in option B.")

# Execute the function to print the results
print_solution()