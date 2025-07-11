import math

def solve_archery_puzzle():
    """
    Solves the archery puzzle to find the number of possible values for the bull's eye score.
    """

    # This helper function checks if a target sum can be achieved by adding k scores from a given list.
    def is_achievable(target, k, scores):
        """
        Checks if `target` can be formed by summing `k` numbers from the `scores` list (with replacement).
        Uses dynamic programming: dp[i][j] is True if sum j can be made with i arrows.
        """
        dp = [[False] * (target + 1) for _ in range(k + 1)]
        dp[0][0] = True

        for i in range(1, k + 1):  # Number of arrows
            for j in range(1, target + 1):  # Target sum
                for s in scores:
                    if j >= s and dp[i - 1][j - s]:
                        dp[i][j] = True
                        break
        return dp[k][target]

    # v_i are the score values divided by 5. They are distinct positive integers.
    # Anna: 3*v1 + v2 + v4 = 25
    # Cliff's simplified score: 185 / 5 = 37
    # Bobby's simplified score: 230 / 5 = 46

    # Constraints derived from player scores: v1 <= 7 and v4 >= 10.
    
    possible_v4_values = set()

    # Iterate through possible values for v1 (from 1 to 7)
    for v1 in range(1, 8):
        sum_v2_v4 = 25 - 3 * v1

        # From v1 < v2 < v3 < v4, we know v4 >= v2 + 2.
        # So, sum_v2_v4 - v2 >= v2 + 2, which simplifies to v2 <= (sum_v2_v4 - 2) / 2.
        v2_upper_bound = math.floor((sum_v2_v4 - 2) / 2)

        # Iterate through possible values for v2
        for v2 in range(v1 + 1, v2_upper_bound + 1):
            v4 = sum_v2_v4 - v2

            # Ensure v4 meets the constraint v4 >= 10
            if v4 < 10:
                continue

            # Check if there exists a valid v3 for this (v1, v2, v4) combination
            is_scenario_valid = False
            for v3 in range(v2 + 1, v4):
                current_scores = [v1, v2, v3, v4]
                
                # Check if Cliff's score (37) is achievable
                if not is_achievable(37, 5, current_scores):
                    continue

                # Check if Bobby's score (46) is achievable
                if not is_achievable(46, 5, current_scores):
                    continue

                # If both scores are achievable, this v4 is possible.
                is_scenario_valid = True
                break  # A valid v3 was found, no need to check others for this (v1, v2, v4)

            if is_scenario_valid:
                possible_v4_values.add(v4)
    
    # The final result is the number of unique possible values for v4.
    # The original bull's eye score is s4 = 5 * v4.
    # The number of possible v4 values is the same as the number of possible s4 values.
    print(f"The simplified possible values for the bull's eye score (divided by 5) are: {sorted(list(possible_v4_values))}")
    print(f"The number of possible values for the score of the bull's eye is: {len(possible_v4_values)}")

solve_archery_puzzle()