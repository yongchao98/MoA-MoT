import sys

def identify_molecule():
    """
    Identifies and provides details about the molecule shown in the image.
    """
    description = "The image shows a large, synthetic macrocycle with C2 rotational symmetry."
    
    # Constituent units
    naphthalene_units = 2
    phenylene_units = 6 # (2 para, 4 meta)
    ethenylene_units = 2
    ethynylene_units = 6
    
    # Calculate chemical formula
    carbon_count = (naphthalene_units * 10) + (phenylene_units * 6) + (ethenylene_units * 2) + (ethynylene_units * 2)
    hydrogen_count = (naphthalene_units * 6) + (phenylene_units * 4) + (ethenylene_units * 2)

    descriptive_name = "Shape-persistent arylene-ethynylene-ethenylene macrocycle"
    chemical_formula = f"C{carbon_count}H{hydrogen_count}"
    origin = "This molecule is 'template 4' from Isobe et al., Nature Chemistry 9, 899–905 (2017)."

    print("Molecule Identification:")
    print(f"Description: {description}")
    print(f"Descriptive Name: {descriptive_name}")
    print(f"Chemical Formula: {chemical_formula}")
    print("\n--- Formula Calculation ---")
    print(f"Number of Carbon atoms = {carbon_count}")
    print(f"Number of Hydrogen atoms = {hydrogen_count}")
    print(f"-------------------------\n")
    print(f"Scientific Reference: {origin}")

if __name__ == "__main__":
    identify_molecule()