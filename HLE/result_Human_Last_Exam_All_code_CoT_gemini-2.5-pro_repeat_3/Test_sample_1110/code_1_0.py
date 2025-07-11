def analyze_boronic_acid_synthesis(start_material_moles, buLi_eq, borate_eq, procedure_name):
    """
    Analyzes the stoichiometry of a boronic acid synthesis to calculate
    the ratio of product to byproduct.

    Args:
        start_material_moles (float): Moles of the starting aryl halide.
        buLi_eq (float): Equivalents of n-BuLi used.
        borate_eq (float): Equivalents of trimethyl borate used.
        procedure_name (str): Name of the procedure being analyzed (e.g., "Original").
    """
    print(f"--- Analysis of the {procedure_name} Procedure ---")

    # The amount of product formed is limited by the starting material, assuming full conversion.
    product_moles = start_material_moles
    
    # Total moles of trimethyl borate added to the reaction
    borate_total_moles = start_material_moles * borate_eq
    
    # Moles of trimethyl borate that react to form the desired product
    borate_reacted_moles = start_material_moles
    
    # Moles of excess trimethyl borate remaining after the reaction
    borate_excess_moles = borate_total_moles - borate_reacted_moles
    
    # During aqueous workup, the excess borate hydrolyzes to form boric acid byproduct
    boric_acid_byproduct_moles = borate_excess_moles
    
    print(f"Assuming we start with {start_material_moles:.2f} mole of 2-bromo-4-chloro-1-iodobenzene.")
    print(f"The final mixture after workup contains:")
    
    # The "final equation" is the composition of the product mixture.
    # We output each number in this final composition.
    print(f"  - Desired Product (Ar-B(OH)2): {product_moles:.2f} moles")
    print(f"  - Boric Acid Byproduct (B(OH)3): {boric_acid_byproduct_moles:.2f} moles")
    
    print("\nResulting Issue:")
    if boric_acid_byproduct_moles > product_moles:
        print("The mixture contains a large excess of boric acid byproduct relative to the desired product.")
        print("This leads to two distinct signals in the Boron NMR and makes purification difficult.")
    elif boric_acid_byproduct_moles > 0:
        print("The amount of boric acid byproduct is significantly reduced.")
        print("This will result in a much purer product, likely showing one major signal in the Boron NMR.")
    else:
        print("No boric acid byproduct is formed, which is ideal.")

# --- Main Script ---

# Let's assume a basis of 1.0 mole of starting material for calculation.
start_moles = 1.0

# Define the conditions from the problem statement
original_buLi_eq = 1.05
original_borate_eq = 5.0

# Define the proposed improved conditions (Answer D: use less trimethyl borate)
# A 1.2x excess is a more reasonable amount.
improved_borate_eq = 1.2

print("Problem: Two Boron NMR signals were observed in a synthesis.\n")

# Analyze the original procedure
analyze_boronic_acid_synthesis(start_moles, original_buLi_eq, original_borate_eq, "Original")

print("\n" + "="*70 + "\n")

# Analyze the proposed solution
print("Proposed Solution: Use less trimethyl borate (Answer D).\n")
analyze_boronic_acid_synthesis(start_moles, original_buLi_eq, improved_borate_eq, "Improved")
