import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string stream
sys.stdout = captured_output = io.StringIO()

def solve_chemistry_problem():
    """
    This function formalizes the logical steps to identify Compound A based on the reaction and spectral data.
    """
    # Step 1: Deconstruct the product from NMR data.
    # The ¹H NMR spectrum shows key fragments.
    # From reactant benzylamine:
    nh_proton = 8.69  # Triplet, 1H
    ch2_protons = 4.73  # Doublet, 2H
    phenyl_protons = "7.37-7.22" # Multiplet, 5H
    
    # From reactant tert-butyl hydrazine:
    tert_butyl_protons = 1.70  # Singlet, 9H
    
    # Protons from the core of Compound A:
    core_proton_1 = 8.24  # Singlet, 1H
    core_proton_2 = 8.11  # Singlet, 1H

    # Step 2: Deduce the core structure.
    # The presence of two singlet protons (8.24 and 8.11 ppm) on the aromatic core is the key clue.
    # This indicates a ring system where two ring protons have no adjacent proton neighbors.
    # A pyrimidine ring substituted at positions 4 and 6 fits this description perfectly,
    # as the remaining protons at positions 2 and 5 are isolated from each other.
    
    # Step 3: Determine the starting material.
    # The reaction is a double nucleophilic substitution on this pyrimidine core.
    # This requires two good leaving groups at the 4 and 6 positions.
    # Halogens are standard leaving groups for such reactions.
    # The most common starting material for this transformation is 4,6-dichloropyrimidine.
    
    compound_a = "4,6-Dichloropyrimidine"
    
    print("The starting material, Compound A, is identified as:")
    print(compound_a)
    
    # As requested, output the numbers from the problem description that led to the conclusion.
    print("\nThis conclusion is based on key signals from the ¹H NMR spectrum of the final product:")
    print(f"An NH proton as a triplet at {nh_proton} ppm.")
    print(f"Two core aromatic protons as singlets at {core_proton_1} ppm and {core_proton_2} ppm.")
    print(f"Two CH2 protons as a doublet at {ch2_protons} ppm.")
    print(f"Nine tert-butyl protons as a singlet at {tert_butyl_protons} ppm.")

solve_chemistry_problem()

# Restore original stdout and print the captured output
sys.stdout = original_stdout
# We use print to ensure the final output is directly displayed to the user
print(captured_output.getvalue())