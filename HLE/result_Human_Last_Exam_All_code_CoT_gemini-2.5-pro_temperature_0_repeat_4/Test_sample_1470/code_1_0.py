import sys

def propose_methanol_charges():
    """
    This function proposes a reasonable set of partial charges for a methanol molecule,
    verifies its neutrality, and prints the results.
    The chosen set corresponds to option D, which is the most chemically sound and
    adheres to the principle of net neutrality and symmetry.
    """
    # Based on analysis, option D is the most reasonable choice.
    # C: 0.1450, H(methyl): 0.0400, O: -0.6830, H(hydroxyl): 0.4180
    charges = {
        "Carbon": 0.1450,
        "Methyl H1": 0.0400,
        "Methyl H2": 0.0400,
        "Methyl H3": 0.0400,
        "Oxygen": -0.6830,
        "Hydroxyl hydrogen": 0.4180
    }

    print("Proposed reasonable partial charge assignments for methanol:")
    # Using sys.stdout.write for cleaner output control without implicit newlines
    for atom, charge in charges.items():
        # Format for alignment: Atom name left-aligned in 20 chars, charge right-aligned in 8
        sys.stdout.write(f"{atom:<20}: {charge:>8.4f}\n")

    print("\nVerifying that the molecule is electrically neutral:")
    
    total_charge = 0
    equation_parts = []
    
    # Iterate through the defined charges to build the equation string and calculate the sum
    for charge in charges.values():
        total_charge += charge
        # Add parentheses for negative numbers for clarity in the equation
        if charge < 0:
            equation_parts.append(f"({charge:.4f})")
        else:
            equation_parts.append(f"{charge:.4f}")
            
    equation_str = " + ".join(equation_parts)
    
    # Print the full equation showing each charge
    print(f"Sum of charges = {equation_str} = {total_charge:.4f}")

if __name__ == "__main__":
    propose_methanol_charges()