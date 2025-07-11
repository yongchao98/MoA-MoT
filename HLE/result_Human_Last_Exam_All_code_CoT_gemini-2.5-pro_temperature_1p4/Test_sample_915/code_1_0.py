def find_combination(total_score, scores, num_arrows):
    """
    Checks if a total_score can be achieved with num_arrows given the score values.
    Returns the combination of hits if possible, otherwise None.
    """
    s1, s2, s3, s4 = scores
    # n1, n2, n3, n4 are the number of hits on each zone
    for n1 in range(num_arrows + 1):
        for n2 in range(num_arrows - n1 + 1):
            for n3 in range(num_arrows - n1 - n2 + 1):
                n4 = num_arrows - n1 - n2 - n3
                # This condition is implicitly handled by the loop ranges, but good for clarity
                if n4 >= 0:
                    current_score = n1*s1 + n2*s2 + n3*s3 + n4*s4
                    if current_score == total_score:
                        return (n1, n2, n3, n4)
    return None

def solve_archery_puzzle():
    """
    Finds the number of possible values for the bull's eye score.
    """
    possible_s1_values = set()
    bobby_score = 230
    cliff_score = 185
    num_arrows = 5

    # From 3*s4 + s3 + s1 = 125 and s1 > s3 > s4 > 0, we can deduce 5*s4 < 125, so s4 < 25.
    # s4 must be a multiple of 5.
    for s4 in range(5, 25, 5):
        # From 3*s4 + s3 + s1 = 125, we have s3 < 125 - 3*s4.
        # s3 must be a multiple of 5 and greater than s4.
        for s3 in range(s4 + 5, 125 - 3 * s4, 5):
            # Calculate s1 from Anna's score equation.
            s1 = 125 - 3 * s4 - s3

            # Constraint 1: Scores must be strictly increasing. s1 > s3.
            if s1 <= s3:
                continue

            # Constraint 2: Bobby's average score is 230/5 = 46.
            # Max score possible for him is 5*s1. So 5*s1 >= 230 -> s1 >= 46.
            # Since s1 is a multiple of 5, s1 must be >= 50.
            if s1 < 50:
                continue

            # Now, iterate through possible s2 values. s3 < s2 < s1.
            for s2 in range(s3 + 5, s1, 5):
                scores = (s1, s2, s3, s4)

                # Check if this score set works for Bobby
                bobby_hits = find_combination(bobby_score, scores, num_arrows)
                if not bobby_hits:
                    continue

                # Check if this score set works for Cliff
                cliff_hits = find_combination(cliff_score, scores, num_arrows)
                if not cliff_hits:
                    continue
                
                # If we are here, we found a valid set of scores.
                # We print the details only for the first time we find a new s1 value.
                if s1 not in possible_s1_values:
                    print(f"--- Found a possible bull's eye score: s1 = {s1} ---")
                    print(f"This is supported by the score set (s1,s2,s3,s4) = {scores}")
                    print("Verification:")
                    print(f"Anna:  3 * {s4} + 1 * {s3} + 1 * {s1} = {3*s4+s3*1+s1*1}")
                    print(f"Bobby: {bobby_hits[0]}*{s1} + {bobby_hits[1]}*{s2} + {bobby_hits[2]}*{s3} + {bobby_hits[3]}*{s4} = {sum(h*s for h,s in zip(bobby_hits, scores))}")
                    print(f"Cliff: {cliff_hits[0]}*{s1} + {cliff_hits[1]}*{s2} + {cliff_hits[2]}*{s3} + {cliff_hits[3]}*{s4} = {sum(h*s for h,s in zip(cliff_hits, scores))}\n")
                
                possible_s1_values.add(s1)

    print("--- Calculation Complete ---")
    print(f"The possible values for the bull's eye score are: {sorted(list(possible_s1_values))}")
    print(f"The exact number of possible values for the score of the bull's eye is: {len(possible_s1_values)}")

solve_archery_puzzle()