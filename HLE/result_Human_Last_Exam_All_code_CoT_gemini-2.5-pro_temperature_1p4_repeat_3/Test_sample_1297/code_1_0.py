import sys
import io

# A function to generate the output
def generate_molecule_info():
    """
    This function provides the SMILES string and molecular weight calculation for the designed molecule.
    """
    
    # The designed molecule is bis(2-morpholinoethyl) ether.
    smiles_string = "O(CCN1CCOCC1)CCN2CCOCC2"
    
    # The molecular formula derived from the problem's numerical constraints is C12H24N2O3.
    # We use precise isotope masses for the calculation to match the target MW.
    # Carbon (12C): 12.00000 Da
    # Hydrogen (1H): 1.00783 Da
    # Nitrogen (14N): 14.00307 Da
    # Oxygen (16O): 15.99491 Da

    num_C = 12
    mass_C = 12.00000
    num_H = 24
    mass_H = 1.00783
    num_N = 2
    mass_N = 14.00307
    num_O = 3
    mass_O = 15.99491

    # Calculate the total molecular weight
    total_mw = (num_C * mass_C) + (num_H * mass_H) + (num_N * mass_N) + (num_O * mass_O)

    # Print the results
    print(f"SMILES: {smiles_string}")
    print("Molecular Weight Calculation:")
    
    # As requested, output each number in the final equation
    calculation_string = f"{num_C} * {mass_C} + {num_H} * {mass_H} + {num_N} * {mass_N} + {num_O} * {mass_O} = {total_mw:.5f}"
    print(calculation_string)

# Execute the function to print the output
generate_molecule_info()
