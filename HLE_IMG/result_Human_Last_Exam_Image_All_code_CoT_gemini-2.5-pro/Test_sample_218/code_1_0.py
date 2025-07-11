def solve_reaction():
    """
    This function identifies the product of the given two-step reaction.
    """
    # Define starting material and product properties
    starting_material = {
        "name": "Geraniol",
        "iupac_name": "(2E)-3,7-dimethylocta-2,6-dien-1-ol",
        "formula": "C10H18O"
    }

    product_A = {
        "name": "Compound A (Dihydromyrcene)",
        "iupac_name": "(6E)-2,6-dimethylocta-2,6-diene",
        "formula": "C10H18"
    }

    # Print the information
    print("Starting Material:")
    for key, value in starting_material.items():
        print(f"  {key.replace('_', ' ').capitalize()}: {value}")

    print("\nReaction Type: Deoxygenation (R-OH -> R-H)")
    
    print("\nFinal Product:")
    for key, value in product_A.items():
        print(f"  {key.replace('_', ' ').capitalize()}: {value}")

# Execute the function to display the result
solve_reaction()