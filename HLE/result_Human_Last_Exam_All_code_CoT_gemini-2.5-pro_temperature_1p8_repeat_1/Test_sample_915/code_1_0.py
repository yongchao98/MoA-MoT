def solve_archery_challenge():
    """
    This function solves the archery challenge by finding all possible score values for the bull's eye
    and returns the count of these possibilities.
    """

    def check_combination(scores, total_score, num_arrows=5):
        """
        Checks if a total score can be achieved with 5 arrows for a given set of ring scores.
        It iterates through all partitions of 5 arrows among the 4 score regions.
        """
        s1, s2, s3, s4 = scores
        for n4 in range(num_arrows + 1):  # Hits on bull's eye
            for n3 in range(num_arrows - n4 + 1):  # Hits on inner ring
                for n2 in range(num_arrows - n4 - n3 + 1):  # Hits on next ring
                    n1 = num_arrows - n4 - n3 - n2  # Hits on outer ring
                    
                    calculated_score = n1 * s1 + n2 * s2 + n3 * s3 + n4 * s4
                    if calculated_score == total_score:
                        return True
        return False

    possible_s4_values = set()

    # Step 1: Iterate through possible scores based on Anna's data and constraints.
    # Anna's equation: 3*s1 + s2 + s4 = 125.
    # Scores s1, s2, s3, s4 are strictly increasing multiples of 5.
    # From 3*s1 < 125, we know s1 < 41.6.
    for s1 in range(5, 45, 5):
        
        # We need space for s3 between s2 and s4, so s4 >= s2 + 10.
        # Substituting into Anna's equation: s2 + (s2 + 10) <= 125 - 3*s1
        # This gives an upper bound for s2: 2*s2 <= 115 - 3*s1.
        s2_upper_bound = (115 - 3 * s1) // 2
        for s2 in range(s1 + 5, s2_upper_bound + 1, 5):
            
            # s4 is determined by Anna's equation.
            s4 = 125 - 3 * s1 - s2
            
            # Iterate through all possible s3 values.
            for s3 in range(s2 + 5, s4, 5):
                scores = (s1, s2, s3, s4)
                
                # Step 2: Check if the score set is valid for Bobby (230).
                if check_combination(scores, 230):
                    
                    # Step 3: Check if the score set is valid for Cliff (185).
                    if check_combination(scores, 185):
                        # If the score set works for everyone, s4 is a possible value.
                        possible_s4_values.add(s4)

    # Step 4: Output the results.
    final_count = len(possible_s4_values)
    sorted_values = sorted(list(possible_s4_values))

    print(f"Anna's equation: 3 * s1 + s2 + s4 = 125")
    print(f"Bobby's score: 230")
    print(f"Cliff's score: 185")
    print("-" * 30)
    print(f"Found {final_count} possible value(s) for the bull's eye score.")
    print(f"The possible values are: {sorted_values}")

solve_archery_challenge()
<<<8>>>