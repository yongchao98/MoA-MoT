def count_homology_cobordism_elements():
    """
    Calculates the number of homology cobordism group elements representable
    by integral surgery on knots with at most four crossings.
    """
    # Define knots with <= 4 crossings. Note: -3_1 is the mirror of 3_1.
    # 4_1 is amphichiral (its own mirror).
    knots = {
        "0_1": {"name": "Unknot", "type": "slice"},
        "3_1": {"name": "Right-handed Trefoil", "type": "chiral"},
        "-3_1": {"name": "Left-handed Trefoil", "type": "chiral_mirror"},
        "4_1": {"name": "Figure-eight knot", "type": "slice"},
    }

    # Surgery on a knot K yields a homology sphere if the framing is +/- 1.
    surgeries = [-1, 1]

    # Use a set to store the unique element identifiers found.
    unique_elements = set()
    element_descriptions = {}

    print("Analyzing integral surgeries on knots with at most four crossings:")
    print("-" * 65)

    # Trivial Element
    trivial_id = "Trivial Element (Class of S^3)"
    print("1. The Trivial Element:")
    print("   - Surgery on the Unknot (0_1) or any slice knot (like 4_1) with +/-1 framing")
    print("     results in a homology sphere that is homology cobordant to S^3.")
    print("     This represents the trivial element in the homology cobordism group.")
    unique_elements.add(trivial_id)
    element_descriptions[trivial_id] = 1 # for final calculation

    # Poincaré Sphere
    poincare_id = "Poincaré Sphere (Class of Sigma(2,3,5))"
    print("\n2. The Poincaré Sphere Element:")
    print("   - Resulting from (-1)-surgery on the Right-handed Trefoil (3_1).")
    print("   - This element has order 2, meaning it is its own inverse.")
    print("   - (+1)-surgery on the Left-handed Trefoil (-3_1) yields the same element.")
    unique_elements.add(poincare_id)
    element_descriptions[poincare_id] = 1 # for final calculation

    # Sigma(2,3,7) and its inverse
    sigma_id = "Brieskorn Sphere (Class of Sigma(2,3,7))"
    sigma_inverse_id = "Inverse of Sigma(2,3,7)"
    print("\n3. The Sigma(2,3,7) Element and its Inverse:")
    print(f"   - (+1)-surgery on the Right-handed Trefoil (3_1) yields the '{sigma_id}'.")
    print("     This element has infinite order.")
    print(f"   - (-1)-surgery on the Left-handed Trefoil (-3_1) yields the '{sigma_inverse_id}'.")
    unique_elements.add(sigma_id)
    unique_elements.add(sigma_inverse_id)
    # This pair contributes 2 to the final sum
    element_descriptions[sigma_id] = 2 # Representing the pair

    print("\nSummary of unique elements found:")
    for i, elem in enumerate(sorted(list(unique_elements))):
        print(f"  - {elem}")

    count = len(unique_elements)
    
    # Build and print the final equation
    calculation_parts = [str(v) for k, v in sorted(element_descriptions.items())]
    calculation_str = " + ".join(calculation_parts)
    
    print("\nFinal Calculation:")
    print(f"The total number of distinct elements is the sum of:")
    print(f"- 1 (for the trivial element)")
    print(f"- 1 (for the Poincaré sphere class)")
    print(f"- 2 (for the Sigma(2,3,7) class and its distinct inverse)")
    print("\nEquation:")
    print(f"{calculation_str} = {count}")

if __name__ == "__main__":
    count_homology_cobordism_elements()