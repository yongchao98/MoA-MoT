def solve_nanoribbon_classification():
    """
    Analyzes and classifies three graphene nanoribbon band structures.
    """
    
    # Classification data derived from visual analysis of the plots
    classifications = [
        {'id': 1, 'Edge': 'A', 'Width': 8, 'Band': 1},
        {'id': 2, 'Edge': 'A', 'Width': 5, 'Band': 0},
        {'id': 3, 'Edge': 'Z', 'Width': 7, 'Band': 0},
    ]

    final_string = ""
    print("Classification Analysis:")
    for item in classifications:
        edge_code = item['Edge']
        width_val = item['Width']
        band_code = item['Band']
        
        # Build the classification string for the current plot
        classification_str = f"{edge_code}{width_val}{band_code}"
        
        # Explain the components for clarity
        edge_text = "Armchair" if edge_code == 'A' else "Zigzag"
        band_text = "semiconducting" if band_code == 1 else "metallic"
        
        print(f"\nPlot {item['id']}:")
        print(f"  - Edge: {edge_code} ({edge_text})")
        print(f"  - Width (N): {width_val}")
        print(f"  - Band: {band_code} ({band_text})")
        print(f"  - Resulting code: {classification_str}")

        # Append to the final concatenated string
        final_string += classification_str

    print("\n" + "="*30)
    print("Final Concatenated Classification:")
    print(final_string)
    print("="*30)

solve_nanoribbon_classification()