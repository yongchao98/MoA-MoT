import sympy

def calculate_conductance():
    """
    This function symbolically calculates the four-terminal conductance G_12,34
    for a quantum Hall device with a QPC, as described in the problem.
    """
    # Define symbolic variables for the number of channels.
    # M is the total number of edge states.
    # N is the number of reflected edge states.
    M, N = sympy.symbols('M N', positive=True, integer=True)

    # Define symbols for elementary charge (e) and Planck's constant (h).
    e, h = sympy.symbols('e h', positive=True, real=True)

    # The quantum of conductance is G_0 = e^2/h.
    G_0 = e**2 / h

    # Based on the Landauer-Büttiker formalism, the relationship between the
    # injected current I (from terminal 1 to 2) and the source-drain
    # voltage (V1 - V2) is derived as:
    # I = G_0 * (M - N) * (V1 - V2)

    # The relationship for the measured voltage difference (V3 - V4) is:
    # V34 = (N / M) * (V1 - V2)

    # The four-terminal conductance G_12,34 is the ratio I / V34.
    # G_12,34 = [G_0 * (M - N) * (V1 - V2)] / [(N / M) * (V1 - V2)]
    # The (V1 - V2) terms cancel out.
    
    conductance_coefficient = (M * (M - N)) / N
    G_12_34 = conductance_coefficient * G_0

    # --- Output ---
    print("This script calculates the symbolic formula for the four-terminal conductance G_12,34.")
    print("The system is a 6-terminal Hall bar with M edge states and a QPC reflecting N states.")
    print("\nVariable Definitions:")
    print(f"  M: Total number of spin-degenerate edge states.")
    print(f"  N: Number of edge states reflected by the QPC.")
    print(f"  e: Elementary charge.")
    print(f"  h: Planck's constant.")
    print(f"  G_0 = e^2/h is the quantum of conductance.")

    print("\nThe final formula for the conductance G_12,34 is derived using the Landauer-Büttiker formalism.")
    print("The result is:")
    
    # Print the equation with the symbolic components M and N, as requested.
    print("\n  G_12,34 = (M * (M - N) / N) * (e^2/h)")
    
    # For a prettier mathematical representation:
    from sympy.printing.pretty import pretty
    pretty_G = pretty(G_12_34, use_unicode=False)
    print("\nSymbolic expression from sympy:")
    print(f"  G_12,34 = {pretty_G}")


if __name__ == '__main__':
    calculate_conductance()
