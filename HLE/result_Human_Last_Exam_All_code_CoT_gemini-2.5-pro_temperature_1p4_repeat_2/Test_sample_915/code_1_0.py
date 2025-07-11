def solve_archery_puzzle():
    """
    Solves the archery puzzle to find the number of possible values for the bull's eye score.
    """

    # A cache to store results of is_formable to avoid re-computation.
    memo = {}

    def is_formable(total, scores, num_shots):
        """
        Checks if a 'total' can be formed by summing 'num_shots' from the 'scores' list.
        Uses dynamic programming.
        """
        # Create a hashable key for the memoization cache
        scores_tuple = tuple(sorted(scores))
        if (total, scores_tuple, num_shots) in memo:
            return memo[(total, scores_tuple, num_shots)]

        # dp[i][j] is True if a sum of j can be made with i shots.
        dp = [[False] * (total + 1) for _ in range(num_shots + 1)]
        dp[0][0] = True

        for i in range(1, num_shots + 1):
            for j in range(1, total + 1):
                for score in scores_tuple:
                    if j >= score:
                        if dp[i - 1][j - score]:
                            dp[i][j] = True
                            break
        
        result = dp[num_shots][total]
        memo[(total, scores_tuple, num_shots)] = result
        return result

    # Set to store the unique possible values for k4 (the bull's eye score scaled by 5).
    possible_k4_values = set()

    # Anna's equation: 3*k1 + k2 + k4 = 25.
    # Constraints: 0 < k1 < k2 < k3 < k4.
    # This implies k2 >= k1+1 and k4 >= k2+2.
    # From these, we can deduce 1 <= k1 <= 4.
    for k1 in range(1, 5):
        # Determine the upper bound for k2 based on Anna's equation and inequalities.
        limit_k2 = (25 - 3 * k1 - 2) // 2 
        for k2 in range(k1 + 1, limit_k2 + 1):
            k4 = 25 - 3 * k1 - k2

            # Check if there's room for k3, i.e., k4 > k2 + 1
            if k4 <= k2 + 1:
                continue

            # Bobby's average score is 230/5 = 46. s4 must be >= 46.
            # s4 = 5*k4, so 5*k4 >= 46 => k4 >= 9.2. So, k4 must be at least 10.
            if k4 < 10:
                continue
            
            # Now, for the valid (k1, k2, k4) triplet, check if a valid k3 exists.
            for k3 in range(k2 + 1, k4):
                scores = (k1, k2, k3, k4)
                
                # Check if Cliff's scaled score (185/5 = 37) is possible.
                cliff_possible = is_formable(37, scores, 5)
                
                if cliff_possible:
                    # If Cliff's score is possible, check Bobby's (230/5 = 46).
                    bobby_possible = is_formable(46, scores, 5)
                    
                    if bobby_possible:
                        # If both are possible, this is a valid scenario.
                        # The bull's eye score k4 is a possible value.
                        possible_k4_values.add(k4)
                        # We only need one valid k3, so we can stop searching.
                        break
                        
    # Print the problem's constants and the final result
    print(f"Anna's equation: 3 * (outer ring) + 1 * (next ring) + 1 * (bull's eye) = 125")
    print(f"Bobby's total score: 230")
    print(f"Cliff's total score: 185")
    print(f"Number of arrows shot by each person: 5")
    print("-" * 20)
    print(f"The number of possible values for the bull's eye score is:")
    print(len(possible_k4_values))

solve_archery_puzzle()