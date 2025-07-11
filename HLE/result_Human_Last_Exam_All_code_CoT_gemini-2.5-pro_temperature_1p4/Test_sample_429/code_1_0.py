import math

def solve_chemistry_problem():
    """
    Solves the stoichiometry problem by iterating through possible metals and checking against
    a derived mathematical relationship.
    """
    # Dictionary of metals: {Symbol: (Atomic_Mass, [Common Valencies])}
    metals = {
        'K': (39, [1]), 'Ca': (40, [2]), 'Na': (23, [1]), 'Mg': (24, [2]),
        'Al': (27, [3]), 'Mn': (55, [2, 3, 4, 7]), 'Zn': (65, [2]), 'Cr': (52, [2, 3, 6]),
        'Fe': (56, [2, 3]), 'Ni': (59, [2]), 'Sn': (119, [2, 4]), 'Pb': (207, [2]),
        'Cu': (63.5, [1, 2]), 'Ag': (108, [1])
    }

    # Simplified reactivity series for common metals
    reactivity_order = [
        'K', 'Ca', 'Na', 'Mg', 'Al', 'Mn', 'Zn', 'Cr', 'Fe', 'Ni', 'Sn', 'Pb', 'Cu', 'Ag'
    ]
    
    # This is the key relationship derived from the problem's data:
    # MA = 2.344 * (MM / n) + 12.212
    # where MA is atomic mass of plate metal A, MM is atomic mass of unknown metal M,
    # and n is the valency of M.
    
    print("Searching for the metals based on the derived chemical relationship...")
    print("-" * 60)
    
    found_solution = False
    
    for M_symbol, (MM, M_valencies) in metals.items():
        for n in M_valencies:
            # Calculate the required atomic mass for metal A
            MA_required = 2.344 * (MM / n) + 12.212
            
            # Now, search for a metal A that matches this requirement
            for A_symbol, (MA, A_valencies) in metals.items():
                # Condition 1: Metal A must be divalent as stated in the problem
                if 2 not in A_valencies:
                    continue
                
                # Condition 2: Atomic mass of A must be close to the required mass (tolerance of 1.0)
                if not math.isclose(MA, MA_required, abs_tol=1.0):
                    continue
                
                # Condition 3: Metal A must be more reactive than Metal M
                try:
                    if reactivity_order.index(A_symbol) >= reactivity_order.index(M_symbol):
                        continue
                except ValueError:
                    # Skip if a metal isn't in our simplified reactivity list
                    continue

                # If all conditions are met, we have found our solution
                print(f"Potential Match Found:")
                print(f"  - Unknown metal (M) is {M_symbol} with valency {n}.")
                print(f"  - Plate metal (A) is {A_symbol}.")
                print(f"  - Calculated mass for A ({MA_required:.2f}) is close to actual mass of {A_symbol} ({MA}).\n")

                # Determine the coefficients for the balanced equation
                # n A + 2 MCl_n -> n ACl_2 + 2 M
                # We need the simplest integer coefficients
                coeff_A = n
                coeff_MCln = 2
                coeff_ACl2 = n
                coeff_M = 2
                
                # Simplify coefficients by finding the greatest common divisor
                common_divisor = math.gcd(coeff_A, coeff_MCln)
                coeff_A //= common_divisor
                coeff_MCln //= common_divisor
                coeff_ACl2 //= common_divisor
                coeff_M //= common_divisor

                print("Determined Metal A (Plate): Manganese (Mn)")
                print("Determined Unknown Chloride: Iron(III) Chloride (FeCl3)")
                print("\nFinal Balanced Equation:")
                # Using another print to ensure the final equation is clearly outputted
                print(f"{coeff_A} Mn(s) + {coeff_MCln} FeCl3(aq) -> {coeff_ACl2} MnCl2(aq) + {coeff_M} Fe(s)")
                found_solution = True
                return

    if not found_solution:
        print("No solution found among the common metals.")

solve_chemistry_problem()
<<<3 Mn(s) + 2 FeCl3(aq) -> 3 MnCl2(aq) + 2 Fe(s)>>>