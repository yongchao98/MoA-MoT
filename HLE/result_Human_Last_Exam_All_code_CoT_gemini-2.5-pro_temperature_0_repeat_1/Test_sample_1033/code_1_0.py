def solve():
    """
    This function calculates the next three triplets in the sequence.
    """
    # The last triplet given is NZX.
    # N = 13, Z = 25, X = 23
    v2_prev = 25
    v3_prev = 23

    # The next group of triplets will start with 'O'.
    # O = 14
    v1_current = 14

    # The equations for the sequence are:
    # v2_current = (v2_prev + v3_prev) % 26
    # v3_current = (v1_current + v2_current) % 26

    print("The rule is determined to be:")
    print("L2_current = (L2_previous + L3_previous) mod 26")
    print("L3_current = (L1_current + L2_current) mod 26")
    print("\nStarting from the last known triplet NZX (N=13, Z=25, X=23):")
    
    results = []
    for i in range(3):
        # Calculate the next triplet
        v2_current = (v2_prev + v3_prev) % 26
        v3_current = (v1_current + v2_current) % 26

        # Convert numerical values back to letters
        l1_current_char = chr(ord('A') + v1_current)
        l2_current_char = chr(ord('A') + v2_current)
        l3_current_char = chr(ord('A') + v3_current)
        
        triplet_str = f"{l1_current_char}{l2_current_char}{l3_current_char}"
        results.append(triplet_str)

        print(f"\n--- Triplet {i+1} ---")
        print(f"L1 = {l1_current_char} ({v1_current})")
        print(f"L2 = ({v2_prev} + {v3_prev}) mod 26 = {v2_current} ({l2_current_char})")
        print(f"L3 = ({v1_current} + {v2_current}) mod 26 = {v3_current} ({l3_current_char})")
        print(f"Result: {triplet_str}")

        # Update the previous values for the next iteration
        v2_prev = v2_current
        v3_prev = v3_current

    print("\nThe next three capital letters in the sequence are:")
    print(" ".join(results))

solve()
<<<OWK OGU OAO>>>