import sys

def propose_methanol_charges():
    """
    Proposes and validates a reasonable set of partial charges for an all-atom
    model of methanol (CH3OH) based on common physical and chemical principles.
    """
    
    print("Proposing a new set of partial charges for an all-atom model of methanol (CH3OH).\n")
    print("The proposed charges are based on the following principles:")
    print("1. The molecule must be charge-neutral.")
    print("2. Chemically equivalent atoms (the three methyl hydrogens) should have identical charges.")
    print("3. Charges should reflect electronegativity differences (e.g., O is negative, H on O is positive).")
    print("-" * 30)

    # Based on analysis, option D is the most reasonable choice.
    # It is charge-neutral, respects symmetry, and has a chemically intuitive charge distribution
    # similar to established force fields like OPLS-AA without being an exact copy.
    charges = {
        'Carbon': 0.1450,
        'Methyl H1': 0.0400,
        'Methyl H2': 0.0400,
        'Methyl H3': 0.0400,
        'Oxygen': -0.6830,
        'Hydroxyl hydrogen': 0.4180
    }

    print("Proposed Partial Charges (in fundamental charge units, e):")
    
    total_charge = 0
    for atom, charge in charges.items():
        # Using a format specifier to ensure charges are aligned and displayed with consistent precision
        print(f"{atom:<18}: {charge:>8.4f}")
        total_charge += charge
        
    print("-" * 30)
    # The sum should be exactly zero. We use a small tolerance for floating point comparisons.
    # f-string formatting with +.7f will show the sign and enough decimal places to confirm neutrality.
    print(f"Verification of Net Charge: {total_charge:+.7f}")

    if abs(total_charge) < 1e-9:
        print("The total charge is neutral, as required for a simulation of a neutral molecule.")
    else:
        print("Warning: The molecule is not charge-neutral!", file=sys.stderr)

propose_methanol_charges()
<<<D>>>