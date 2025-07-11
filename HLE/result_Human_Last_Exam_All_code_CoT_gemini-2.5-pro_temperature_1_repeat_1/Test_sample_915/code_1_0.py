def solve_archery_puzzle():
    """
    Solves the archery puzzle to find the number of possible values for the bull's eye score.
    """

    def is_achievable(scores, total_score):
        """
        Checks if a total_score can be achieved with 5 arrows given the four possible scores.
        Iterates through all combinations of 5 shots.
        """
        s1, s2, s3, sb = scores
        # nb, n3, n2, n1 are the number of hits on the bull's eye, ring 3, ring 2, and ring 1.
        for nb in range(6):  # Number of hits on bull's eye (0 to 5)
            for n3 in range(6 - nb):  # Number of hits on ring 3
                for n2 in range(6 - nb - n3):  # Number of hits on ring 2
                    n1 = 5 - nb - n3 - n2  # Remaining hits are on ring 1
                    
                    current_score = (n1 * s1) + (n2 * s2) + (n3 * s3) + (nb * sb)
                    
                    if current_score == total_score:
                        return True
        return False

    possible_sb_values = set()

    # 1. Iterate through possible values for s1, the outermost ring score.
    # From the plan, we know s1 can only be 5, 10, 15, or 20.
    for s1 in range(5, 25, 5):

        # 2. Iterate through possible values for s2.
        # From Anna's equation: s2 + sb = 125 - 3*s1
        # From the ordering rule s1 < s2 < s3 < sb, we know sb >= s2 + 10.
        # So, s2 + (s2 + 10) <= 125 - 3*s1  =>  2*s2 <= 115 - 3*s1
        s2_upper_bound = (115 - 3 * s1) / 2
        for s2 in range(s1 + 5, int(s2_upper_bound) + 1, 5):

            # 3. Calculate sb from Anna's equation
            sb = 125 - 3 * s1 - s2
            
            # 4. Check if there's space for a valid s3.
            # We need s2 < s3 < sb, so the gap between sb and s2 must be at least 10.
            if (sb - s2) >= 10:
                
                # 5. Iterate through all possible values for s3.
                for s3 in range(s2 + 5, sb, 5):
                    scores = (s1, s2, s3, sb)

                    # Optimization: check if scores are in a possible range for Bobby and Cliff
                    bobby_possible = (5 * s1 <= 230) and (5 * sb >= 230)
                    cliff_possible = (5 * s1 <= 185) and (5 * sb >= 185)

                    if not (bobby_possible and cliff_possible):
                        continue

                    # 6. Check if this set of scores works for Bobby (230) and Cliff (185)
                    if is_achievable(scores, 230) and is_achievable(scores, 185):
                        possible_sb_values.add(sb)
    
    # Final Output
    print(f"Anna's score can be represented by the equation: 3*s1 + s2 + sb = 125")
    print(f"Applying all constraints, the possible scores for the bull's eye are: {sorted(list(possible_sb_values))}")
    print(f"The exact number of possible values for the score of the bull's eye is: {len(possible_sb_values)}")

solve_archery_puzzle()
<<<3>>>