import math

def check_sum_possible(target_sum, num_shots, scores):
    """
    Checks if a target_sum can be achieved by picking num_shots from the scores list.
    Uses dynamic programming to solve this subset-sum-like problem.
    dp[i][j] will be True if a sum of j is possible using i shots.
    """
    dp = [[False] * (target_sum + 1) for _ in range(num_shots + 1)]
    dp[0][0] = True

    for i in range(1, num_shots + 1):
        for j in range(1, target_sum + 1):
            for score in scores:
                if j >= score and dp[i-1][j-score]:
                    dp[i][j] = True
                    break
    return dp[num_shots][target_sum]

def find_possible_bullseye_scores():
    """
    Finds all possible values for the bull's eye score based on the problem's constraints.
    """
    possible_b_values = set()
    anna_total = 25  # 125 / 5
    bobby_total = 46  # 230 / 5
    cliff_total = 37  # 185 / 5

    # Iterate through possible values for the outer ring score (r3)
    # From 3*r3 + r2 + b = 25 and b > r2 > r3 >= 1, we can deduce 5*r3 < 23, so r3 can be 1, 2, 3, 4.
    for r3 in range(1, 5):
        # Iterate through possible values for the next ring score (r2)
        # From b > r2, we get 25 - r2 - 3*r3 > r2 => r2 < (25 - 3*r3)/2
        r2_upper_bound = (anna_total - 3 * r3) / 2
        for r2 in range(r3 + 1, int(r2_upper_bound) + 1):
            # The upper bound is strict, so we check again
            if r2 >= r2_upper_bound:
                continue

            # Calculate the bull's eye score (b) from Anna's equation
            b = anna_total - r2 - 3 * r3

            # The scores must be strictly increasing, so b > r2.
            # We also need space for r1, so b > r1 > r2, which means b must be at least r2 + 2.
            if b <= r2 + 1:
                continue

            # A potential b is found. Now check if there is at least one valid r1.
            is_b_possible = False
            for r1 in range(r2 + 1, b):
                scores = [b, r1, r2, r3]

                # Condition 1: The GCD of all scores must be 1.
                if math.gcd(math.gcd(b, r1), math.gcd(r2, r3)) != 1:
                    continue

                # Condition 2: Check if Bobby's score (46) is possible with 5 shots.
                if not check_sum_possible(bobby_total, 5, scores):
                    continue

                # Condition 3: Check if Cliff's score (37) is possible with 5 shots.
                if not check_sum_possible(cliff_total, 5, scores):
                    continue

                # If all conditions are met, this b is a possible value.
                # We can stop checking other r1 values for this b.
                is_b_possible = True
                break
            
            if is_b_possible:
                possible_b_values.add(b)

    # Convert the simplified 'b' values back to the actual score values by multiplying by 5
    final_scores = sorted([val * 5 for val in possible_b_values])
    
    print(f"Anna's score equation: 3 * s_r3 + 1 * s_r2 + 1 * s_b = 125")
    print(f"Bobby's score with 5 arrows: 230")
    print(f"Cliff's score with 5 arrows: 185")
    print("-" * 30)
    print(f"The possible values for the bull's eye score are: {final_scores}")
    print(f"The exact number of possible values for the score of the bull's eye is: {len(final_scores)}")

if __name__ == "__main__":
    find_possible_bullseye_scores()