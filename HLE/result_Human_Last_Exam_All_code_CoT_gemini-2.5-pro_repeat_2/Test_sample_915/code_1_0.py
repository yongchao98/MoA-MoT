import collections

# memoization table for the recursive checker
# The key is a tuple (target_score, num_arrows, scores_tuple)
# to ensure the memo is specific to the score set being checked.
memo = {}

def is_score_possible(target_score, num_arrows, scores):
    """
    Checks if a target_score can be achieved by summing num_arrows scores
    chosen from the 'scores' list. Uses recursion with memoization.
    """
    # Normalize scores by dividing by 5 to work with smaller numbers.
    target_norm = target_score // 5
    scores_norm = tuple(sorted([s // 5 for s in scores]))

    # Create a unique key for the memoization table.
    memo_key = (target_norm, num_arrows, scores_norm)
    if memo_key in memo:
        return memo[memo_key]

    # Base cases for the recursion.
    if target_norm == 0 and num_arrows == 0:
        return True
    if target_norm < 0 or num_arrows == 0:
        return False
        
    # Recursive step: try subtracting each possible score.
    for s_norm in scores_norm:
        if is_score_possible(target_score - (s_norm * 5), num_arrows - 1, scores):
            memo[memo_key] = True
            return True
    
    memo[memo_key] = False
    return False

def find_possible_bullseye_scores():
    """
    Main function to find the number of possible values for the bull's eye score.
    """
    possible_bullseye_scores = set()

    # Anna's equation: 3*x1 + x2 + x4 = 25, where x_i = s_i / 5.
    # We know 0 < x1 < x2 < x4.
    # From x1 < x2, we get 3*x1 + x1 < 3*x1 + x2 => 4*x1 < 25-x4.
    # The smallest possible x1 is 1. 5*x1 < 25 => x1 < 5.
    for x1 in range(1, 5):
        # From x2 < x4, we get x2 < 25 - 3*x1 - x2 => 2*x2 < 25 - 3*x1.
        upper_bound_x2 = (25 - 3 * x1) / 2
        for x2 in range(x1 + 1, int(upper_bound_x2)):
            x4 = 25 - 3 * x1 - x2

            # The constraint s2 < s3 < s4 implies x2 < x3 < x4.
            # This means there must be an integer between x2 and x4, so x4 - x2 > 1.
            if x4 - x2 <= 1:
                continue

            s1, s2, s4 = 5 * x1, 5 * x2, 5 * x4
            
            # Check if there exists an s3 that validates this (s1, s2, s4) triplet.
            # A triplet is valid if there's an s3 for which both Bobby's and Cliff's scores are possible.
            found_valid_s3 = False
            for x3 in range(x2 + 1, x4):
                s3 = 5 * x3
                current_scores = [s1, s2, s3, s4]
                
                # Clear memo for each new set of scores being tested
                memo.clear()
                bobby_ok = is_score_possible(230, 5, current_scores)
                
                if bobby_ok:
                    memo.clear()
                    cliff_ok = is_score_possible(185, 5, current_scores)
                    if cliff_ok:
                        found_valid_s3 = True
                        break  # A valid s3 was found, no need to check others for this triplet.
            
            if found_valid_s3:
                possible_bullseye_scores.add(s4)

    return len(possible_bullseye_scores)

if __name__ == '__main__':
    # Calculate and print the result.
    num_possible_values = find_possible_bullseye_scores()
    print(num_possible_values)
    # The problem has a specific format for the final answer.
    print(f"<<<{num_possible_values}>>>")
