from collections import Counter

def format_formula(atoms, name):
    """Formats a dictionary of atom counts into a chemical formula string."""
    # Standard order of elements in a formula: C, H, then alphabetical for the rest.
    formula = ""
    if 'C' in atoms and atoms['C'] > 0:
        formula += f"C{atoms['C']}"
    if 'H' in atoms and atoms['H'] > 0:
        formula += f"H{atoms['H']}"
    
    other_elements = sorted([el for el in atoms if el not in ['C', 'H']])
    
    for el in other_elements:
        if atoms[el] > 0:
            formula += f"{el}{atoms[el]}"
            
    print(f"{name}: {formula}")

def main():
    """
    Calculates the molecular formula of the product through the reaction sequence.
    """
    # Step 0: Define the starting material's molecular formula.
    # The bicyclic core is C7H5F3NO. The PMB protecting group is C8H9O.
    # Starting Material = Core + PMB group = C(7+8)H(5+9)F3N1O(1+1) = C15H14F3NO2
    start_material = Counter({'C': 15, 'H': 14, 'F': 3, 'N': 1, 'O': 2})
    format_formula(start_material, "Starting Material")

    # Step 1: Deprotection of PMB group with CAN.
    # This reaction removes the PMB group (C8H9O) and adds one Hydrogen (H).
    print("\nStep 1: CAN, ACN/H2O (Deprotection of PMB)")
    pmb_group = Counter({'C': 8, 'H': 9, 'O': 1})
    hydrogen_atom = Counter({'H': 1})
    intermediate_1 = start_material - pmb_group + hydrogen_atom
    format_formula(intermediate_1, "Intermediate 1")

    # Step 2: Hydrogenation with Pd/C, H2.
    # This reaction reduces the C=C double bond, adding one molecule of H2.
    print("\nStep 2: Pd/C, H2, MeOH (Hydrogenation)")
    h2_molecule = Counter({'H': 2})
    intermediate_2 = intermediate_1 + h2_molecule
    format_formula(intermediate_2, "Intermediate 2")

    # Step 3: Acid hydrolysis of the lactam.
    # This reaction opens the lactam ring by adding one molecule of H2O.
    print("\nStep 3: 4 N HCl, 70 C (Lactam Hydrolysis)")
    water_molecule = Counter({'H': 2, 'O': 1})
    final_product = intermediate_2 + water_molecule
    format_formula(final_product, "Final Product")
    
    print("\nThe molecular formula of the final product is C7H10F3NO2.")
    print("C=7, H=10, F=3, N=1, O=2")

if __name__ == "__main__":
    main()
