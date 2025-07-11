import sys
from itertools import combinations_with_replacement

def solve_archery_puzzle():
    """
    This script finds the number of possible values for the bull's eye score based on the problem's constraints.
    It systematically searches for all valid sets of scores and verifies them against the information provided for all three archers.
    """
    valid_bull_eye_scores = set()

    # Scores are sB, s1, s2, s3. They are multiples of 5.
    # Let b, r1, r2, r3 be the integer multipliers (score = 5 * multiplier).
    # Rule 1: b > r1 > r2 > r3 >= 1
    # Anna's score equation: 3*s3 + s2 + sB = 125  =>  3*r3 + r2 + b = 25

    # Iterate through possible values for the score multipliers, starting from the outermost ring.
    # The range for r3 is limited. Since b > r1 > r2 > r3, at a minimum, r2=r3+1, r1=r3+2, b=r3+3.
    # Plugging into the equation: 3*r3 + (r3+1) + (r3+3) <= 25 => 5*r3 <= 21 => r3 <= 4.2
    for r3 in range(1, 5):  # r3 can be 1, 2, 3, 4.
        
        # The range for r2 is also constrained.
        # We need to find an integer r1 such that b > r1 > r2, so b must be at least r2 + 2.
        # b = 25 - 3*r3 - r2.
        # So, 25 - 3*r3 - r2 > r2 + 1 => 24 - 3*r3 > 2*r2 => r2 < (24 - 3*r3) / 2.
        r2_upper_bound = int((24 - 3 * r3 - 1) / 2)
        for r2 in range(r3 + 1, r2_upper_bound + 1):
            
            # The bull's eye multiplier 'b' is now determined by r3 and r2.
            b = 25 - 3 * r3 - r2
            sB = b * 5
            
            # Now iterate through possible r1 values, which must be between r2 and b.
            for r1 in range(r2 + 1, b):
                # A complete, potential set of scores is now defined.
                s3 = r3 * 5
                s2 = r2 * 5
                s1 = r1 * 5
                scores = [sB, s1, s2, s3]

                # Check if this score set works for Bobby (230) and Cliff (185).
                # Generate all possible total scores for 5 arrows.
                possible_sums = {sum(combo) for combo in combinations_with_replacement(scores, 5)}
                
                # If both totals are achievable, this is a valid set of scores.
                if 230 in possible_sums and 185 in possible_sums:
                    valid_bull_eye_scores.add(sB)
                    # Once a valid r1 is found for a given (r3,r2,b), sB is confirmed as
                    # a possibility. We can break to avoid redundant checks with other r1 values.
                    break
    
    # The final answer is the number of unique possible values found for the bull's eye score.
    final_count = len(valid_bull_eye_scores)
    print(final_count)

solve_archery_puzzle()