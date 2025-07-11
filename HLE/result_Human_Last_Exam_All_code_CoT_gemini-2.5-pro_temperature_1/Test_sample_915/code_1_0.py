def is_achievable(target_score, scores, num_arrows=5):
    """
    Checks if a target_score can be achieved by summing num_arrows from the given scores.
    The scores tuple is (s1, s2, s3, s_bull).
    This function uses nested loops to check all combinations of 5 arrows.
    """
    s1, s2, s3, s_bull = scores
    
    # An optimization: check if the score is even possible.
    min_possible_score = num_arrows * s1
    max_possible_score = num_arrows * s_bull
    if not (min_possible_score <= target_score <= max_possible_score):
        return False

    for n1 in range(num_arrows + 1):
        for n2 in range(num_arrows - n1 + 1):
            for n3 in range(num_arrows - n1 - n2 + 1):
                n4 = num_arrows - n1 - n2 - n3
                # The loop ranges ensure n4 is always 0 or positive.
                
                current_score = n1 * s1 + n2 * s2 + n3 * s3 + n4 * s_bull
                if current_score == target_score:
                    return True
    return False

def solve_archery_puzzle():
    """
    Finds the number of possible values for the bull's eye score.
    """
    possible_bull_eye_scores = set()

    # From Anna's score, we have 3a + b + c = 25, where s1=5a, s2=5b, s_bull=5c.
    # Constraints: 1 <= a < b and there must be a score s3 between s2 and s_bull.
    # The s3 constraint implies b+1 < c.

    # Iterate through possible values for a (s1/5). 'a' can only be 1, 2, 3, or 4.
    for a in range(1, 5):
        # Iterate through possible values for b (s2/5).
        for b in range(a + 1, 25):
            c = 25 - 3 * a - b
            
            # Apply the constraints: b must be less than c, and there must be
            # room for an intermediate score s3.
            if b + 1 < c:
                s1 = 5 * a
                s2 = 5 * b
                s_bull = 5 * c
                
                # This (s1, s2, s_bull) is a candidate. Now we check if there exists
                # an s3 that makes Bobby's and Cliff's scores possible.
                found_valid_s3 = False
                for s3_candidate in range(s2 + 5, s_bull, 5):
                    scores = (s1, s2, s3_candidate, s_bull)
                    
                    bobby_is_possible = is_achievable(230, scores)
                    if bobby_is_possible:
                        cliff_is_possible = is_achievable(185, scores)
                        if cliff_is_possible:
                            # We found a valid s3, so this s_bull is possible.
                            found_valid_s3 = True
                            break  # Exit the s3 loop, no need to check other s3 values.
                
                if found_valid_s3:
                    possible_bull_eye_scores.add(s_bull)

    # Print the final result per the user instructions.
    final_count = len(possible_bull_eye_scores)
    print(f"The possible values for the bull's eye are: {sorted(list(possible_bull_eye_scores))}")
    print(f"The final equation is: Number of possible values = {final_count}")

solve_archery_puzzle()