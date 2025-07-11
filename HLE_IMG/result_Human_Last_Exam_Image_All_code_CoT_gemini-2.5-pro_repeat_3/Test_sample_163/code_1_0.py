def calculate_molecular_weight(formula):
    """Calculates the molecular weight of a compound from its chemical formula."""
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'O': 15.999,
    }
    
    import re
    # Find all element-count pairs in the formula
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    
    total_weight = 0
    for element, count in elements:
        count = int(count) if count else 1
        if element in atomic_weights:
            total_weight += atomic_weights[element] * count
        else:
            print(f"Warning: Atomic weight for element {element} not found.")
            return None
            
    return total_weight

def main():
    """
    Identifies the products of the reaction and prints their details.
    """
    product_A_name = "2-(tert-butoxy)-1-phenylethyl benzoate"
    product_B_name = "1-(tert-butoxy)-2-phenylethyl benzoate"
    
    # Both products are isomers and have the same molecular formula.
    # Styrene (C8H8) + tert-butyl peroxybenzoate (C11H14O3) -> Product (C19H22O3)
    molecular_formula = "C19H22O3"
    
    molecular_weight = calculate_molecular_weight(molecular_formula)
    
    print("The reaction produces two major products, A and B, which are constitutional isomers.")
    print("-" * 60)
    
    print("Product A:")
    print(f"  Name: {product_A_name}")
    
    print("\nProduct B:")
    print(f"  Name: {product_B_name}")
    
    print("-" * 60)
    print("Properties for both isomers (A and B):")
    print(f"  Molecular Formula: {molecular_formula}")
    print(f"  Molecular Weight: {molecular_weight:.3f} g/mol")

    # As requested, printing each number from the final "equation" (molecular formula)
    print("\nNumbers from the molecular formula (C, H, O):")
    c_count = 19
    h_count = 22
    o_count = 3
    print(f"  Number of Carbon atoms: {c_count}")
    print(f"  Number of Hydrogen atoms: {h_count}")
    print(f"  Number of Oxygen atoms: {o_count}")


if __name__ == '__main__':
    main()
