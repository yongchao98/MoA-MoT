def design_liquid_crystal(n):
    """
    Generates and prints the details for a liquid crystal molecule
    based on the alkyl chain length 'n'.
    
    The general structure is C(n)H(2n+1)-Ph-CN.
    """
    
    # Calculate the number of hydrogen atoms in the alkyl chain
    hydrogen_count = 2 * n + 1
    
    # Generate the chemical formula string
    formula = f"C{n}H{hydrogen_count}-Ph-CN"
    
    print(f"Alkyl Chain Length (n): {n}")
    print(f"Resulting Molecular Formula: {formula}")
    print("Numbers in the final equation (C(n)H(2*n+1)):")
    # As requested, output each number in the equation part CnH2n+1
    print(f"n = {n}")
    print(f"2 = 2")
    print(f"n = {n}")
    print(f"1 = 1")
    print("-" * 30)

def main():
    """
    Main function to run the liquid crystal design simulation.
    """
    print("Liquid Crystal Design Assistant")
    print("Goal: Design a single-ring molecule with room temperature transitions.")
    print("General Structure: C(n)H(2n+1)-Ph-CN\n")
    
    # Step 1: Start with the recommended base molecule (n=5)
    print("Step 1: Initial design recommendation (n=5)")
    print("This molecule, 4-pentylbenzonitrile, is known to have a nematic phase near room temperature.")
    design_liquid_crystal(5)
    
    # Step 2: Demonstrate tuning for a high transition temperature
    print("Step 2: Tuning Scenario - If transition temperature is TOO HIGH")
    print("Strategy: Increase the alkyl chain length to increase molecular flexibility.")
    print("Suggestion: Try n=6.")
    design_liquid_crystal(6)
    
    # Step 3: Demonstrate tuning for a low transition temperature
    print("Step 3: Tuning Scenario - If transition temperature is TOO LOW")
    print("Strategy: Decrease the alkyl chain length to decrease flexibility.")
    print("Suggestion: Try n=4.")
    design_liquid_crystal(4)

if __name__ == "__main__":
    main()