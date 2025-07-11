def solve_archery_puzzle():
    """
    Solves the archery puzzle by systematically checking all possible score combinations
    based on the given constraints and reports the number of possible values for the
    bull's eye score.
    """

    # Memoization dictionary to store results for check_sums_possible,
    # avoiding re-computation for the same set of scores.
    memo_check = {}

    def check_sums_possible(scores_tuple):
        """
        Checks if the scores for Cliff (185) and Bobby (230) can be formed
        by summing 5 scores from the given tuple using dynamic programming.
        """
        scores_tuple = tuple(sorted(scores_tuple))
        if scores_tuple in memo_check:
            return memo_check[scores_tuple]

        # reachable_sums[k] will store the set of all possible sums using k shots.
        reachable_sums = {0: {0}}
        
        for k in range(1, 6):
            current_sums = set()
            for score in scores_tuple:
                for prev_sum in reachable_sums[k - 1]:
                    current_sums.add(prev_sum + score)
            reachable_sums[k] = current_sums

        cliff_possible = 185 in reachable_sums[5]
        bobby_possible = 230 in reachable_sums[5]
        
        result = cliff_possible and bobby_possible
        memo_check[scores_tuple] = result
        return result

    valid_bullseye_scores = set()

    # Iterate through all possible values for the bull's eye score (s4).
    # From Bobby's score, 5*s4 >= 230 => s4 >= 46, so s4 starts at 50.
    # From Anna's score, s4 is at most 100.
    for s4 in range(50, 101, 5):
        is_s4_value_possible = False
        
        # Iterate through possible values for the outer ring score (s1).
        # From Cliff's score, 5*s1 <= 185 => s1 <= 37, so s1 is at most 35.
        for s1 in range(5, 36, 5):
            # Calculate s2 from Anna's equation: 3*s1 + s2 + s4 = 125.
            s2 = 125 - (3 * s1) - s4

            # s2 must be a multiple of 5 and strictly greater than s1.
            if s2 % 5 != 0 or s2 <= s1:
                continue

            # Iterate through possible values for s3.
            for s3 in range(s2 + 5, s4, 5):
                scores = (s1, s2, s3, s4)

                # Check if this score set allows for Bobby's and Cliff's totals.
                if check_sums_possible(scores):
                    valid_bullseye_scores.add(s4)
                    is_s4_value_possible = True
                    # Break as we only need one valid scenario for this s4.
                    break 
            
            if is_s4_value_possible:
                break
    
    # The final question asks for the exact number of possible values.
    # Let's print the possible values and their count.
    # To satisfy the "final equation" instruction, we show the numbers that add up to the count.
    # The sum is figurative in this case.
    if valid_bullseye_scores:
        equation_str = " + ".join(["1" for _ in valid_bullseye_scores])
        print(f"The possible values for the bull's eye are: {sorted(list(valid_bullseye_scores))}")
        print(f"The number of possibilities is given by the sum: {equation_str} = {len(valid_bullseye_scores)}")
    else:
        print("No possible values for the bull's eye score were found.")
        print("The number of possibilities is: 0")
    
# Run the solver
solve_archery_puzzle()