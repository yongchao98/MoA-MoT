def solve_rope_puzzle():
    """
    Solves the rope puzzle by simulating the cutting process.
    """
    initial_length = 60
    shrink_factor = 0.75
    target_length = 15
    solution_found = False

    # 1. Iterate through all possible first cuts.
    # The longer piece `len1` must be > 30cm.
    for len1 in range(initial_length // 2 + 1, initial_length):
        if solution_found:
            break
        
        len2 = initial_length - len1

        # The longer piece `len1` shrinks. New length must be an integer.
        shrunk_len1 = len1 * shrink_factor
        if shrunk_len1 != int(shrunk_len1):
            continue
        shrunk_len1 = int(shrunk_len1)

        # After the first cut, the two available pieces are `shrunk_len1` and `len2`.
        # We can perform the second cut on either of these.
        pieces_for_cut2 = [shrunk_len1, len2]
        
        for rope_for_cut2 in pieces_for_cut2:
            if solution_found:
                break
            
            # 2. Iterate through all possible second cuts on the chosen piece.
            for cut2_len1 in range(1, rope_for_cut2):
                cut2_len2 = rope_for_cut2 - cut2_len1
                
                # Identify the longer and shorter pieces of the second cut.
                if cut2_len1 > cut2_len2:
                    longer2, shorter2 = cut2_len1, cut2_len2
                else:
                    longer2, shorter2 = cut2_len2, cut2_len1

                # The longer piece shrinks. New length must be an integer.
                shrunk_len2 = longer2 * shrink_factor
                if shrunk_len2 != int(shrunk_len2):
                    continue
                shrunk_len2 = int(shrunk_len2)
                
                # 3. Check if we achieved the target length.
                # The final piece can't be the one that shrunk, because if it was,
                # its original length would have been 15 / 0.75 = 20, making the other piece smaller,
                # which contradicts that the 20cm piece was the longer one.
                # So, the target piece must be `shorter2`.
                if shorter2 == target_length:
                    # We found the solution.
                    print("Solution found through the following steps:")
                    print(f"1. First, the 60cm rope is cut into {len1}cm and {len2}cm pieces.")
                    print(f"   The longer piece ({len1}cm) shrinks to {shrunk_len1}cm.")
                    print(f"2. Second, the resulting {rope_for_cut2}cm piece is cut.")
                    print("The final equation to get the 15cm rope is:")
                    # Output the equation for the second cut
                    print(f"{shorter2} = {rope_for_cut2} - {longer2}")
                    
                    # Determine the longer portion after the first cut.
                    longer_portion_after_cut1 = max(shrunk_len1, len2)
                    print("\nThe length of the longer portion after the first cut is:")
                    print(longer_portion_after_cut1)
                    
                    solution_found = True
                    break

# Run the solver
solve_rope_puzzle()