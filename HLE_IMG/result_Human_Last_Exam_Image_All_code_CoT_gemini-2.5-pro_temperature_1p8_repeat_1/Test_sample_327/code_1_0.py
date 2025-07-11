def solve_graphene_classification():
    """
    Solves the classification problem for the three graphene nanoribbon band structures.
    """
    
    # Classification based on visual analysis of the band structure plots.
    # ribbon = {'id': id, 'Edge': 'A' or 'Z', 'Width': N, 'Band': 0 or 1}
    
    # Plot 1: Armchair (gap at k=0), Semiconducting (gap exists), Width=7 (7 bands above E=0, consistent with semiconducting AGNR).
    ribbon1_class = {'id': 1, 'Edge': 'A', 'Width': 7, 'Band': 1}
    
    # Plot 2: Armchair (linear crossing at k=0), Metallic (no gap), Width=5 (5 bands above E=0, consistent with metallic AGNR).
    ribbon2_class = {'id': 2, 'Edge': 'A', 'Width': 5, 'Band': 0}
    
    # Plot 3: Zigzag (flat band at E=0), Metallic (no gap), Width=6 (6 bands above E=0).
    ribbon3_class = {'id': 3, 'Edge': 'Z', 'Width': 6, 'Band': 0}
    
    ribbons = [ribbon1_class, ribbon2_class, ribbon3_class]
    
    final_string = ""
    
    for ribbon in ribbons:
        classification_str = f"{ribbon['Edge']}{ribbon['Width']}{ribbon['Band']}"
        print(f"Ribbon {ribbon['id']}: Edge={ribbon['Edge']}, Width={ribbon['Width']}, Band={ribbon['Band']} -> {classification_str}")
        final_string += classification_str
        
    print("\nFinal concatenated classification string:")
    print(final_string)

solve_graphene_classification()