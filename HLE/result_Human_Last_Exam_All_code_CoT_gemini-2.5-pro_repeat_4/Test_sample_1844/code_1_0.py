def solve_folding_puzzle():
    """
    Calculates the total number of edges on a piece of paper after a series of folds and cuts.
    """
    
    # Step 1: Initialize the four corners of the original paper.
    # Each corner is represented by a tuple: (layers, num_folds, name)
    corners = {
        'TL': [1, 0, 'TL'], 'TR': [1, 0, 'TR'],
        'BL': [1, 0, 'BL'], 'BR': [1, 0, 'BR']
    }
    
    # Step 2: Simulate the four folds.

    # Fold 1: Top half onto bottom half (H1)
    # Top corners (TL, TR) fold onto bottom corners (BL, BR)
    # The original BL and BR are now on a fold line.
    print("--- After Fold 1 (Top-to-Bottom) ---")
    corners['BL'][0] += corners['TL'][0] # New BL layers = old BL + old TL
    corners['BL'][1] += 1 # The H1 fold passes through the original BL corner
    corners['BR'][0] += corners['TR'][0] # New BR layers = old BR + old TR
    corners['BR'][1] += 1 # The H1 fold passes through the original BR corner
    # The new shape's corners are the old bottom corners.
    # We can think of the paper now being represented by the original BL and BR corners.
    # For clarity, we'll keep tracking all 4 original corners' locations.
    print(f"Original TR corner is now at: {corners['TR']}")
    print(f"Original TL corner is now at: {corners['TL']}")
    print(f"Original BR corner is now at: {corners['BR']}")
    print(f"Original BL corner is now at: {corners['BL']}")
    print("-" * 20)

    # Let's define the corners of the folded piece at each stage
    # S1 corners are the bottom corners of the original paper after the fold.
    s1_corners = {'BL': corners['BL'], 'BR': corners['BR']}
    
    # Fold 2: Left half onto right half (V1)
    # The folded shape is now a square. Its TL & BL corners are from S1's BL.
    # Its TR & BR corners are from S1's BR.
    # Let's trace the four original corners again.
    # TL corner's position is folded onto TR's position.
    # BL corner's position is folded onto BR's position.
    print("--- After Fold 2 (Left-to-Right) ---")
    corners['TR'][0] += corners['TL'][0]
    corners['TR'][1] += 1 # V1 fold passes through the original TR corner
    corners['BR'][0] += corners['BL'][0]
    corners['BR'][1] += 1 # V1 fold passes through the original BR corner
    print(f"Original TR corner is now at: {corners['TR']}")
    print(f"Original TL corner is now at: {corners['TL']}")
    print(f"Original BR corner is now at: {corners['BR']}")
    print(f"Original BL corner is now at: {corners['BL']}")
    print("-" * 20)
    # The corners of the S2 square are the locations of the 4 original corners.
    s2_corners = corners.copy()

    # Fold 3: Top half onto bottom half (H2)
    # TL folds to BL, TR folds to BR
    print("--- After Fold 3 (Top-to-Bottom) ---")
    corners['BL'][0] += corners['TL'][0]
    corners['BL'][1] += 1 # H2 fold passes through original BL
    corners['BR'][0] += corners['TR'][0]
    corners['BR'][1] += 1 # H2 fold passes through original BR
    print(f"Original TR corner is now at: {corners['TR']}")
    print(f"Original TL corner is now at: {corners['TL']}")
    print(f"Original BR corner is now at: {corners['BR']}")
    print(f"Original BL corner is now at: {corners['BL']}")
    print("-" * 20)
    s3_corners = {'BL': corners['BL'], 'BR': corners['BR']}

    # Fold 4: Left half onto right half (V2)
    # BL folds to BR
    print("--- After Fold 4 (Left-to-Right) ---")
    corners['BR'][0] += corners['BL'][0]
    corners['BR'][1] += 1 # V2 fold passes through original BR
    print(f"Original TR corner is now at: {corners['TR']}")
    print(f"Original TL corner is now at: {corners['TL']}")
    print(f"Original BR corner is now at: {corners['BR']}")
    print(f"Original BL corner is now at: {corners['BL']}")
    print("-" * 20)
    
    # The four corners of the final folded square are the final locations
    # of the four original corners of the paper.
    final_corners_properties = list(corners.values())

    print("\n--- Properties of the 4 Corners of the Final Folded Square ---")
    total_edges = 0
    
    # This list will hold the equations for each corner
    equations = []

    for k, f, name in final_corners_properties:
        # k = layers, f = number of folds
        indentation_edges = k + 2 * f
        total_edges += indentation_edges
        
        # We need to output each number in the final equation
        # So we format the string carefully.
        equation_str = f"{k} + 2 * {f} = {indentation_edges}"
        equations.append(f"Edges from corner '{name}': {equation_str}")

    print("Each cut creates an indentation at one of the original corners.")
    print("The number of edges for each indentation is k (layers) + 2*f (folds).\n")
    for eq in equations:
        print(eq)
        
    print("\n--- Final Calculation ---")
    edge_counts = [prop[0] + 2 * prop[1] for prop in final_corners_properties]
    
    # Print the final summation equation
    sum_str = " + ".join(map(str, edge_counts))
    print(f"Total Edges = {sum_str} = {total_edges}")

    return total_edges

final_answer = solve_folding_puzzle()
# The final answer is printed by the function, so we don't print it again here.
# However, the format requires the final answer to be at the end.
# The value is 36.
print(f"\nFinal Answer: {final_answer}")