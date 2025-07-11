def solve_resistor_field():
    """
    This function prints the derived expressions for the electric field
    in the two regions of the cylindrical half-shell resistor.
    """

    # Define symbolic variables as strings for printing
    V0 = "V_0"
    sigma1 = "σ_1"
    sigma2 = "σ_2"
    r = "r"
    pi = "π"
    i_phi = "î_φ"

    # Expression for the electric field in Region 1 (0 < phi < pi/2)
    # E_1 = (2 * σ_2 * V_0) / (r * π * (σ_1 + σ_2)) * î_φ
    print("The electric field in Region 1 is:")
    e1_numerator = f"2 * {sigma2} * {V0}"
    e1_denominator = f"{r} * {pi} * ({sigma1} + {sigma2})"
    print(f"E_1 = ({e1_numerator}) / ({e1_denominator}) * {i_phi}\n")

    # Expression for the electric field in Region 2 (pi/2 < phi < pi)
    # E_2 = (2 * σ_1 * V_0) / (r * π * (σ_1 + σ_2)) * î_φ
    print("The electric field in Region 2 is:")
    e2_numerator = f"2 * {sigma1} * {V0}"
    e2_denominator = f"{r} * {pi} * ({sigma1} + {sigma2})"
    print(f"E_2 = ({e2_numerator}) / ({e2_denominator}) * {i_phi}\n")

    print("This corresponds to answer choice C.")

if __name__ == "__main__":
    solve_resistor_field()