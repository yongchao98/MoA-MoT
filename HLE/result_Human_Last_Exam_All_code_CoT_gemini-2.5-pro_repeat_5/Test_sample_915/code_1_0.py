import itertools

def solve_archery_puzzle():
    """
    This function calculates the number of possible values for the bull's eye score
    based on the given constraints and scores of the three participants.
    """
    # A set to store the unique, valid bull's eye scores.
    valid_bulls_eye_scores = set()

    # The problem can be simplified by dividing all scores by 5.
    # Let k_i = s_i / 5. Anna's equation becomes: k_b + k_r2 + 3*k_r3 = 25.
    # The score order is k_b > k_r1 > k_r2 > k_r3 >= 1.

    # 1. Iterate through all possible score combinations that satisfy Anna's score.
    # The maximum possible value for k_r3 is 4, otherwise k_b > k_r2 > k_r3 is impossible.
    for k_r3 in range(1, 5):
        s_r3 = k_r3 * 5
        
        # Determine the range for k_r2.
        # k_b + k_r2 = 25 - 3*k_r3. Let S = 25 - 3*k_r3.
        # We need k_b > k_r2, so S - k_r2 > k_r2, which means S > 2*k_r2 or k_r2 < S/2.
        # Also, k_r2 must be greater than k_r3.
        s_val = 25 - 3 * k_r3
        for k_r2 in range(k_r3 + 1, (s_val + 1) // 2):
            s_r2 = k_r2 * 5
            
            # Calculate k_b from the equation.
            k_b = 25 - 3 * k_r3 - k_r2
            s_b = k_b * 5

            # 2. For each valid set (s_b, s_r2, s_r3), check if there exists an s_r1
            #    that makes Bobby's and Cliff's scores possible.
            found_valid_s_r1 = False
            # Iterate through all possible values for the missing score s_r1.
            for s_r1 in range(s_r2 + 5, s_b, 5):
                scores = [s_b, s_r1, s_r2, s_r3]
                
                # 3. Generate all possible sums for 5 arrows.
                # Use itertools.combinations_with_replacement for an efficient way to do this.
                possible_sums = {sum(c) for c in itertools.combinations_with_replacement(scores, 5)}
                
                # 4. Check if both Bobby's and Cliff's scores are achievable.
                if 230 in possible_sums and 185 in possible_sums:
                    found_valid_s_r1 = True
                    break  # Found a working s_r1, this s_b is possible.

            if found_valid_s_r1:
                valid_bulls_eye_scores.add(s_b)

    # 5. The final answer is the number of unique possible values for the bull's eye.
    print(f"The number of possible values for the bull's eye score is: {len(valid_bulls_eye_scores)}")

if __name__ == '__main__':
    solve_archery_puzzle()