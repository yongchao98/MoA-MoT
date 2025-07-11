from itertools import combinations_with_replacement

def is_achievable(target_sum, num_arrows, items):
    """
    Checks if a target sum can be achieved by summing num_arrows 
    from the given items, with replacement.
    """
    for combo in combinations_with_replacement(items, num_arrows):
        if sum(combo) == target_sum:
            return True
    return False

def find_bulls_eye_values():
    """
    Calculates the number of possible values for the bull's eye score
    based on the given archery puzzle rules.
    """
    # All scores are multiples of 5, so we can divide by 5 to simplify.
    anna_total_div5 = 25  # 125 / 5
    bobby_total_div5 = 46  # 230 / 5
    cliff_total_div5 = 37  # 185 / 5
    num_arrows = 5
    
    # Store unique possible values for the bull's eye score (s4).
    possible_s4_scores = set()
    
    # We iterate through possible values for 'a' (outer ring score / 5).
    # From 3a + b + d = 25 and a < b < b+1 < d, we can deduce 5a < 21, so a <= 4.
    for a in range(1, 5):
        s1 = 5 * a
        # Iterate through 'b' (next ring score / 5), where b > a.
        # The upper bound is derived from 3a + b + d = 25 and d >= b+2.
        for b in range(a + 1, 11): 
            s2 = 5 * b
            # Calculate 'd' (bull's eye score / 5) from Anna's equation.
            d = anna_total_div5 - 3 * a - b
            s4 = 5 * d

            # Constraint: score values must be strictly increasing.
            # a < b < c < d requires that d > b + 1. If not, break the inner loop.
            if d <= b + 1:
                break

            # Constraint: Bobby's score is 230. The max possible score (5*s4) must be >= 230.
            # 5 * s4 >= 230  => s4 >= 46 => d >= 9.2 => d must be at least 10.
            if d < 10:
                continue

            # We have a valid triplet (s1, s2, s4). Now check if there exists
            # a valid s3 that makes Bobby's and Cliff's scores possible.
            # Iterate through possible 'c' values, where b < c < d.
            for c in range(b + 1, d):
                s3 = 5 * c
                
                scores_div5 = [a, b, c, d]
                
                # Check if this score set works for Bobby.
                bobby_ok = is_achievable(bobby_total_div5, num_arrows, scores_div5)
                if not bobby_ok:
                    continue # Try the next possible s3
                
                # Check if this score set works for Cliff.
                cliff_ok = is_achievable(cliff_total_div5, num_arrows, scores_div5)
                
                # If both are possible, this s4 is a valid solution.
                if cliff_ok:
                    # To fulfill the "output each number in the final equation" rule,
                    # we show the set of scores that validates this s4 value.
                    print(f"Found valid bull's eye score: {s4}. "
                          f"Possible set of scores (s1,s2,s3,s4): {{{s1}, {s2}, {s3}, {s4}}}")
                    possible_s4_scores.add(s4)
                    # We found a working set, no need to check other s3 values for this s4.
                    break
                    
    print("\n--- Summary ---")
    print(f"All possible values for the bull's eye score are: {sorted(list(possible_s4_scores))}")
    print(f"\nThe exact number of possible values for the score of the bull's eye is: {len(possible_s4_scores)}")

if __name__ == '__main__':
    find_bulls_eye_values()