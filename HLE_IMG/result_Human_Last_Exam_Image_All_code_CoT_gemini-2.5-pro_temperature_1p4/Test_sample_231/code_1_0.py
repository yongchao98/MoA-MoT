import collections

def solve_chemistry_problem():
    """
    This script identifies compound C and calculates its molecular weight.
    """

    # Step 1: Deduce the structure of Compound C based on the reaction scheme.
    # The reaction sequence leads to a demethylated, amino-substituted dibenzo[b,e]pyrylium salt.
    # Structure of C: 1,6,8-trihydroxy-3-diethylamino-10-(2,4,6-trihydroxyphenyl)dibenzo[b,e]pyrylium
    # The molecular formula of this cation is [C23H22NO7]+.

    compound_name = "1,6,8-trihydroxy-3-diethylamino-10-(2,4,6-trihydroxyphenyl)dibenzo[b,e]pyrylium"
    molecular_formula = "[C23H22NO7]+"

    print(f"The final product, compound C, is {compound_name}.")
    print(f"Its molecular formula as a cation is {molecular_formula}.")
    print("\nTo provide a quantitative answer, we will calculate its molecular weight.")

    # Step 2: Define atomic weights and atom counts.
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999
    }
    
    atom_counts = {
        'C': 23,
        'H': 22,
        'N': 1,
        'O': 7
    }

    # Step 3: Calculate the molecular weight, showing the equation with numbers.
    print("\nThe final equation for the molecular weight (MW) is:")
    print("MW = (atomic_weight_C * count_C) + (atomic_weight_H * count_H) + (atomic_weight_N * count_N) + (atomic_weight_O * count_O)")
    
    # Printing each number in the final equation as requested.
    print(f"MW = ({atomic_weights['C']} * {atom_counts['C']}) + ({atomic_weights['H']} * {atom_counts['H']}) + ({atomic_weights['N']} * {atom_counts['N']}) + ({atomic_weights['O']} * {atom_counts['O']})")

    # Perform the calculation
    mw_c = atomic_weights['C'] * atom_counts['C']
    mw_h = atomic_weights['H'] * atom_counts['H']
    mw_n = atomic_weights['N'] * atom_counts['N']
    mw_o = atomic_weights['O'] * atom_counts['O']
    total_mw = mw_c + mw_h + mw_n + mw_o

    print(f"MW = {mw_c:.3f} + {mw_h:.3f} + {mw_n:.3f} + {mw_o:.3f}")
    print(f"\nThe calculated molecular weight of the cation is: {total_mw:.3f} g/mol")

solve_chemistry_problem()