def solve_magnetization_curve():
    """
    This function explains the derivation and prints the analytical expression
    for the initial magnetization curve of a superconducting bar based on the
    Bean critical-state model.
    """

    # Define symbolic variables for clarity in the explanation
    a = "a"   # half-width of the bar
    Jc = "Jc" # critical current density
    H = "H"   # applied magnetic field

    # --- Print the explanation of the model and derivation ---
    print("This problem is solved using the Bean critical-state model for a superconducting slab.")
    print("Based on the geometry (a bar with width 2a >> thickness 2b), the problem is simplified to a 1D model of an infinite slab of thickness 2*a, neglecting demagnetization effects.\n")
    print("Key Steps:")
    print("1. The applied field H is screened by critical currents ±Jc that penetrate a distance 'd' from the surfaces.")
    print("2. The penetration depth 'd' is proportional to the applied field: d = H / Jc.")
    print("3. The magnetization M is defined as the difference between the spatially-averaged internal field <H_int> and the applied field H.")
    print("4. The average internal field is calculated to be: <H_int> = H^2 / (2 * a * Jc).\n")

    # --- Print the final analytical expression ---
    print("The resulting analytical expression for the initial magnetization curve M(H) is:")
    print("-" * 70)
    # The final equation is printed below, showing all numerical coefficients as requested.
    print(f"M(H) = H^2 / (2 * {a} * {Jc}) - 1 * H")
    print("-" * 70)
    print("\nWhere:")
    print("  M(H) is the magnetization.")
    print("  H    is the applied magnetic field (valid for 0 <= H <= a*Jc).")
    print(f"  {a}    is half the width of the bar in the x-direction.")
    print(f"  {Jc}   is the constant critical current density.")

# Execute the function to get the answer
solve_magnetization_curve()