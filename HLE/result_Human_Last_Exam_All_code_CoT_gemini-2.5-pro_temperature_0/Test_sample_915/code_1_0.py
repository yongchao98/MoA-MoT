def solve_archery_puzzle():
    """
    This function solves the archery puzzle by iterating through all
    possible score combinations and checking them against the given conditions.
    """

    # This helper function checks if a total_score can be achieved with num_shots
    # given the scores of the four target areas.
    def can_make_score(scores, total_score, num_shots):
        s1, s2, s3, s4 = scores
        # Iterate through all possible combinations of num_shots shots.
        # n1, n2, n3, n4 are the number of hits on each ring.
        for n4 in range(num_shots + 1):
            for n3 in range(num_shots - n4 + 1):
                for n2 in range(num_shots - n4 - n3 + 1):
                    n1 = num_shots - n4 - n3 - n2
                    if n1 >= 0:
                        current_score = n1 * s1 + n2 * s2 + n3 * s3 + n4 * s4
                        if current_score == total_score:
                            return True  # Found a valid combination
        return False  # No combination found

    possible_bull_eye_scores = set()

    # Let the scores be s1, s2, s3, s4 for the outer ring to the bull's eye.
    # We know s1 < s2 < s3 < s4 and they are all positive multiples of 5.
    # Anna's score: 3*s1 + s2 + s4 = 125.

    # From 3*s1 + s2 + s4 = 125 and s1 < s2 < s4, we can deduce that s1 must be
    # a multiple of 5 less than 25. So, s1 can be 5, 10, 15, 20.
    for s1 in range(5, 25, 5):

        # From s2 + s4 = 125 - 3*s1 and s2 < s4, we have 2*s2 < 125 - 3*s1.
        # Also, s1 < s2.
        s2_start = s1 + 5
        s2_end = (125 - 3 * s1) // 2
        for s2 in range(s2_start, s2_end + 1, 5):

            # Calculate s4 based on Anna's score.
            s4 = 125 - 3 * s1 - s2

            # Now, we must find an s3 such that s2 < s3 < s4.
            # s3 must also be a multiple of 5.
            for s3 in range(s2 + 5, s4, 5):
                scores = (s1, s2, s3, s4)

                # Check if this set of scores is possible for Bobby (total 230)
                bobby_possible = can_make_score(scores, 230, 5)

                if bobby_possible:
                    # If Bobby's score is possible, check Cliff's (total 185)
                    cliff_possible = can_make_score(scores, 185, 5)

                    if cliff_possible:
                        # If both are possible, this is a valid set of scores.
                        # Add the bull's eye score to our set of possibilities.
                        possible_bull_eye_scores.add(s4)

    # The final answer is the number of unique possible values for the bull's eye.
    num_possible_values = len(possible_bull_eye_scores)
    
    print(f"The set of possible scores for the bull's eye is: {sorted(list(possible_bull_eye_scores))}")
    print(f"The number of possible values for the bull's eye is: {num_possible_values}")

solve_archery_puzzle()