import sympy

def calculate_qhe_conductance():
    """
    Symbolically calculates the four-terminal conductance G_12,34 for a
    six-terminal quantum Hall device with a quantum point contact (QPC).
    """
    # Define symbolic variables for the number of channels and fundamental constants
    M, N = sympy.symbols('M N', integer=True, positive=True, doc="Number of channels")
    e, h = sympy.symbols('e h', real=True, positive=True, doc="Fundamental constants")
    V = sympy.Symbol('V', real=True, doc="Source potential")

    # The conductance quantum G0 is defined as e^2/h
    G0 = e**2 / h

    print("--- Derivation of Four-Terminal Conductance G_12,34 ---")
    print(f"Total number of edge states: M")
    print(f"Number of reflected edge states by QPC: N")
    print(f"The conductance is calculated in units of the conductance quantum G0 = e^2/h.")

    # Step 1: Determine the potentials at the voltage probes.
    # We set V_1 = V and V_2 = 0.
    # The potentials of the floating probes are the average of incoming potentials.
    # Edge propagation is clockwise: 1 -> 5 -> 6 -> 2 -> 4 -> 3 -> 1.

    # Potential at terminal 5 (receives M channels from terminal 1 at potential V)
    V_5 = V
    print(f"\nPotential at terminal 5, V_5 = V_1 = V")

    # Potential at terminal 4 (receives M channels from terminal 2 at potential 0)
    V_4 = 0
    print(f"Potential at terminal 4, V_4 = V_2 = 0")

    # Potential at terminal 3 (receives N channels from V_5 and M-N from V_4)
    V_3 = (N * V_5 + (M - N) * V_4) / M
    print(f"Potential at terminal 3, V_3 = (N*V_5 + (M-N)*V_4) / M = {sympy.simplify(V_3)}")

    # Step 2: Calculate the measured voltage difference V_34
    V_34 = V_3 - V_4
    print(f"\nMeasured voltage V_34 = V_3 - V_4 = {sympy.simplify(V_34)}")

    # Step 3: Calculate the source-drain current I
    # Using the Landauer-Büttiker formula for the current at terminal 1:
    # I_1 = G0 * (Sum of incoming currents - Sum of outgoing currents)
    # The channels arriving at terminal 1 come from terminal 3.
    # I_1 = G0 * (M * V_3 - M * V_1)
    I = G0 * (M * V_3 - M * V)
    print(f"Source-drain current I = G0 * (M*V_3 - M*V_1) = {sympy.simplify(I)}")

    # Step 4: Calculate the conductance G_12,34 = I / V_34
    G_12_34 = sympy.simplify(I / V_34)

    # Step 5: Display the final result
    print("\n--- Final Result ---")
    print("The four-terminal conductance G_12,34 is I / V_34.")

    # For clearer output, substitute G0 with a symbol
    G0_sym = sympy.Symbol('G0')
    final_expression = G_12_34.subs(G0, G0_sym)
    print(f"\nG_12,34 = {final_expression}")

    # Break down the final equation as requested
    factor = final_expression / G0_sym
    numerator, denominator = factor.as_numer_denom()

    print("\nBreaking down the final equation:")
    print(f"The conductance is proportional to the conductance quantum G0 = e^2/h.")
    print(f"The proportionality factor is a fraction with:")
    print(f"  - Numerator: {numerator}")
    print(f"  - Denominator: {denominator}")

    print(f"\nSo, the final equation is: G_12,34 = ({numerator}) / ({denominator}) * (e^2/h)")


if __name__ == '__main__':
    calculate_qhe_conductance()