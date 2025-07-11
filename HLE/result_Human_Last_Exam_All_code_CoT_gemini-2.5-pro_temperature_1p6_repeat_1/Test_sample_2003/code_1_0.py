import math

def get_final_formula():
    """
    Derives and simplifies the formula for the sum of sharps.

    This function codifies the step-by-step music theory logic to arrive
    at the final formula for n > 0.
    """
    
    # Step 1: Define the base data for the n=0 chromatic scale.
    # Note Name, Number of initial sharps (k_0), Base sharps for the diatonic note name.
    # Base Sharps (C=0, D=2, E=4, F=-1, G=1, A=3, B=5)
    initial_notes = [
        {'name': 'C',  'k0': 0, 'base_sharps': 0},
        {'name': 'C#', 'k0': 1, 'base_sharps': 0},
        {'name': 'D',  'k0': 0, 'base_sharps': 2},
        {'name': 'D#', 'k0': 1, 'base_sharps': 2},
        {'name': 'E',  'k0': 0, 'base_sharps': 4},
        {'name': 'F',  'k0': 0, 'base_sharps': -1},
        {'name': 'F#', 'k0': 1, 'base_sharps': -1},
        {'name': 'G',  'k0': 0, 'base_sharps': 1},
        {'name': 'G#', 'k0': 1, 'base_sharps': 1},
        {'name': 'A',  'k0': 0, 'base_sharps': 3},
        {'name': 'A#', 'k0': 1, 'base_sharps': 3},
        {'name': 'B',  'k0': 0, 'base_sharps': 5},
    ]

    # Step 2: Calculate the sum of NetSharps for n=0.
    # NetSharps = base_sharps + k0 * 7
    s0_net = sum(note['base_sharps'] + note['k0'] * 7 for note in initial_notes)

    # Step 3: The coefficient for n is 12 * 7, as each of the 12 notes adds 7n sharps.
    n_coefficient = 12 * 7

    # Step 4: The formula for the sum for n > 0 is S(n) = s0_net + n_coefficient * n
    # No adjustment is needed for flats because for n>0, all key signatures have >= 6 sharps.
    
    constant_term = s0_net
    
    # Step 5: Simplify the formula by finding the GCD.
    common_divisor = math.gcd(constant_term, n_coefficient)
    
    simplified_constant = constant_term // common_divisor
    simplified_n_coefficient = n_coefficient // common_divisor
    
    # Output the logic and the final derived formula.
    print("The total sum of sharps S(n) for n > 0 is derived from the formula:")
    print(f"S(n) = Sum_{{i=1}}^{{12}} (BaseSharps_i + (k_{{0,i}} + n) * 7)")
    print(f"S(n) = Sum(BaseSharps_i + k_{{0,i}}*7) + Sum(7n)")
    print(f"S(n) = S_0_net + 12 * 7 * n")
    print("\nCalculating the constant term S_0_net (the sum of net sharps for n=0):")
    print(f"S_0_net = {s0_net}")
    print("\nCalculating the coefficient for n:")
    print(f"Coefficient = 12 notes * 7 sharps per 'n' = {n_coefficient}")
    
    print("\nThus, the formula for n > 0 is:")
    print(f"S(n) = {constant_term} + {n_coefficient}n")
    
    print("\nSimplifying the formula by factoring out the greatest common divisor:")
    print(f"GCD({constant_term}, {n_coefficient}) = {common_divisor}")
    print("Final Simplified Formula for n > 0:")
    print(f"S(n) = {common_divisor} * ({simplified_constant} + {simplified_n_coefficient}n)")

get_final_formula()
# The problem asks for the formula itself as the final answer.
# Based on the derivation, the simplified formula is 6 * (9 + 14n).
# The format should be <<<answer content>>>.
print("\n<<<6 * (9 + 14n)>>>")