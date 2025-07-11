def calculate_product_formula():
    """
    Calculates the molecular formula of the product from the given reaction scheme.
    """
    # Step 1: Determine the molecular formula of the starting material.
    # Core structure: 6-(trifluoromethyl)-2-azabicyclo[2.2.1]hept-5-en-3-one
    # C=7, H=6, F=3, N=1, O=1 (after accounting for H on N)
    # PMB group (p-methoxybenzyl): C8H9O
    # Starting material = Core (with NH) + PMB - H = C7H5F3NO + C8H9O
    c_atoms = 15
    h_atoms = 14
    f_atoms = 3
    n_atoms = 1
    o_atoms = 2
    print("Starting Material Formula: C15H14F3NO2")

    # Step 2: First reaction - Deprotection of PMB group with CAN.
    # The PMB group (C8H9O) is removed and replaced by a Hydrogen atom.
    # Net change is subtraction of C8H8O.
    c_atoms -= 8
    h_atoms -= 8
    o_atoms -= 1
    print(f"After Step 1 (deprotection), Intermediate 1 Formula: C{c_atoms}H{h_atoms}F{f_atoms}N{n_atoms}O{o_atoms}")

    # Step 3: Second reaction - Hydrogenation with Pd/C, H2.
    # The C=C double bond is reduced, adding two Hydrogen atoms.
    h_atoms += 2
    print(f"After Step 2 (hydrogenation), Intermediate 2 Formula: C{c_atoms}H{h_atoms}F{f_atoms}N{n_atoms}O{o_atoms}")

    # Step 4: Third reaction - Hydrolysis of the lactam with HCl.
    # The lactam ring is opened by adding one molecule of water (H2O).
    h_atoms += 2
    o_atoms += 1
    print(f"After Step 3 (hydrolysis), Final Product Formula: C{c_atoms}H{h_atoms}F{f_atoms}N{n_atoms}O{o_atoms}")
    
    # Final Answer
    print("\nThe molecular formula of the product is composed of the following atoms:")
    print(f"Carbon (C): {c_atoms}")
    print(f"Hydrogen (H): {h_atoms}")
    print(f"Fluorine (F): {f_atoms}")
    print(f"Nitrogen (N): {n_atoms}")
    print(f"Oxygen (O): {o_atoms}")
    print("\nFinal Molecular Formula String:")
    print(f"C{c_atoms}H{h_atoms}F{f_atoms}N{n_atoms}O{o_atoms}")


calculate_product_formula()