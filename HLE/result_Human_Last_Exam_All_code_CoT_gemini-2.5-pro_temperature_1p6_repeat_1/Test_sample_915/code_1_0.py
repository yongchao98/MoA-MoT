import collections

def solve():
    """
    Solves the archery puzzle to find the number of possible values for the bull's eye score.
    """

    def get_possible_sums(scores, num_arrows):
        """
        Calculates all possible sums achievable with a given number of arrows and scores.
        """
        # Start with a sum of 0, achievable with 0 arrows.
        possible_sums = {0}
        for _ in range(num_arrows):
            new_sums = set()
            for s in scores:
                for existing_sum in possible_sums:
                    new_sums.add(s + existing_sum)
            possible_sums = new_sums
        return possible_sums

    # --- Main logic ---

    # Set to store the unique valid bull's eye scores (original value, not simplified)
    valid_bull_eye_scores = set()

    # Simplified scores for Bobby and Cliff
    bobby_target_simplified = 230 // 5
    cliff_target_simplified = 185 // 5
    num_arrows = 5

    # Step 1 & 2: Generate candidate tuples (x1, x2, x4) from Anna's data
    # 3*x1 + x2 + x4 = 25
    # Looping x1 from 1. If 3*x1 > 25, no solution is possible.
    for x1 in range(1, 25 // 3):
        # Looping x2 from x1 + 1.
        for x2 in range(x1 + 1, 25):
            x4 = 25 - 3 * x1 - x2
            
            # Since x1 < x2 < x4, if x4 is no longer greater than x2,
            # increasing x2 further will only decrease x4, so we can break.
            if x4 <= x2:
                break
            
            # Constraint: There must be an integer x3 such that x2 < x3 < x4.
            # This implies that x4 must be at least x2 + 2.
            if x4 >= x2 + 2:
                # This is a valid candidate tuple from Anna's info.
                candidate_tuple = (x1, x2, x4)
                
                # Step 3: Iterate through possible x3 values to check for validity.
                is_tuple_valid = False
                for x3 in range(candidate_tuple[1] + 1, candidate_tuple[2]):
                    scores_simplified = [candidate_tuple[0], candidate_tuple[1], x3, candidate_tuple[2]]
                    
                    # Step 4: Check if Bobby's and Cliff's scores are achievable.
                    all_sums = get_possible_sums(scores_simplified, num_arrows)
                    
                    if bobby_target_simplified in all_sums and cliff_target_simplified in all_sums:
                        # Found a valid x3. The bull's eye value x4 is possible.
                        # Store the original score value.
                        bull_eye_score = candidate_tuple[2] * 5
                        valid_bull_eye_scores.add(bull_eye_score)
                        is_tuple_valid = True
                        break # Move to the next candidate tuple
                
    # Step 5: Print the final result.
    num_possible_values = len(valid_bull_eye_scores)
    print(f"The exact number of possible values for the score of the bull's eye is: {num_possible_values}")

solve()