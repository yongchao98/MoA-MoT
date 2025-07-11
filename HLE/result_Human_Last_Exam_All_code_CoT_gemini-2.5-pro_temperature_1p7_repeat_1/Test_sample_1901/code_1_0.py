import sys

def design_liquid_crystal(n):
    """
    This function designs a single benzene ring liquid crystal molecule
    based on a given alkyl chain length 'n'.
    """

    # Key requirements and molecular features for the design
    core_structure = "Ph"  # Represents a single benzene ring
    polar_group = "CN"    # A strong polar group for liquid crystal behavior

    # Print the design plan
    print("--- Liquid Crystal Design Plan ---")
    print("Goal: Design a single benzene ring liquid crystal with a room temperature nematic/smectic phase.")
    print("Strategy: Use the general structure C(n)H(2n+1)-Ph-CN and tune 'n'.\n")

    # Set the alkyl chain length based on the input
    print(f"Step 1: Choose the alkyl chain length 'n'.")
    print(f"Selected value for n = {n}\n")

    # Calculate the number of hydrogen atoms in the alkyl chain
    # The formula is C_n H_{2n+1}
    # Per the instructions, we will output each number in the final equation.
    # The equation is: number_of_hydrogens = 2 * n + 1
    multiplier = 2
    addition = 1
    num_hydrogens = multiplier * n + addition
    
    print(f"Step 2: Calculate the number of hydrogen atoms in the alkyl chain.")
    print(f"The formula for the number of hydrogen atoms is: {multiplier} * n + {addition}")
    print(f"Calculation: {multiplier} * {n} + {addition} = {num_hydrogens}\n")
    
    # Assemble the final molecular structure and print the result
    final_structure = f"C{n}H{num_hydrogens}-{core_structure}-{polar_group}"
    
    print("--- Final Proposed Molecular Structure ---")
    print(f"The general structure is: CnH(2n+1)-Ph-CN")
    print(f"For n={n}, the specific structure is: {final_structure}")
    print("This corresponds to the molecule 4-pentylbenzonitrile.")


# Run the design process for n=5, which is known to be a good starting point for
# room temperature liquid crystals of this type.
design_liquid_crystal(n=5)