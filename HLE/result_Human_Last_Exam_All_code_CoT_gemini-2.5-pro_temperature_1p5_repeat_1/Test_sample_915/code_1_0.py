def solve_archery_puzzle():
    """
    Finds the number of possible values for the bull's eye score based on the
    scores of three archers.
    """
    valid_sb_values = set()

    def can_achieve_total(x_scores, target_x_sum, num_arrows=5):
        """
        Checks if a combination of 'num_arrows' shots with given 'x_scores'
        can result in 'target_x_sum'.
        The scores and sum are divided by 5 to work with smaller integers.
        """
        x1, x2, x3, xb = x_scores
        # Iterate through all possible combinations of num_arrows shots across the 4 score regions
        for n1 in range(num_arrows + 1):
            for n2 in range(num_arrows - n1 + 1):
                for n3 in range(num_arrows - n1 - n2 + 1):
                    nb = num_arrows - n1 - n2 - n3
                    
                    current_x_sum = n1 * x1 + n2 * x2 + n3 * x3 + nb * xb
                    if current_x_sum == target_x_sum:
                        return True
        return False

    # Main logic derived from Anna's score: 3*x1 + x2 + xb = 25
    # The condition 0 < x1 < x2 < x3 < xb is applied.
    
    # We can deduce that x1 must be in {1, 2, 3, 4}.
    for x1 in range(1, 5):
        # We can also deduce a bound for x2 from x1 < x2 and x2 < xb - 1.
        # This leads to 2*x2 < 24 - 3*x1.
        upper_bound_x2 = (24 - 3 * x1) // 2 + 1
        for x2 in range(x1 + 1, upper_bound_x2):
            xb = 25 - 3 * x1 - x2
            
            # We must be able to fit x3 between x2 and xb.
            for x3 in range(x2 + 1, xb):
                
                # A potential set of scores (divided by 5)
                candidate_x_scores = (x1, x2, x3, xb)
                
                # Check if this set works for Bobby (target score 230 -> x_sum 46)
                bobby_possible = can_achieve_total(candidate_x_scores, 46)
                
                if bobby_possible:
                    # If it works for Bobby, check Cliff (target score 185 -> x_sum 37)
                    cliff_possible = can_achieve_total(candidate_x_scores, 37)
                    
                    if cliff_possible:
                        # This score set is valid for all archers.
                        # Add the actual bull's eye score to our set of solutions.
                        bulls_eye_score = 5 * xb
                        valid_sb_values.add(bulls_eye_score)

    # Print the results in a clear format
    sorted_values = sorted(list(valid_sb_values))
    num_values = len(sorted_values)

    print(f"Found {num_values} possible values for the score of the bull's eye.")
    if num_values > 0:
        print(f"The possible values are: {', '.join(map(str, sorted_values))}")
        
        # Per instructions, showing the final "equation" for the count
        equation_str = " + ".join(["1"] * num_values)
        print(f"The final equation for the count of possible values is: {equation_str} = {num_values}")


# Execute the solution
solve_archery_puzzle()