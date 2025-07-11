from itertools import combinations_with_replacement

def solve_archery_challenge():
    """
    This script calculates the number of possible values for the bull's eye score
    based on the information provided for the three archers.
    """

    # Helper function that checks if a total_score can be formed by summing
    # num_shots values from the available_scores list.
    def is_sum_possible(available_scores, total_score, num_shots):
        for combo in combinations_with_replacement(available_scores, num_shots):
            if sum(combo) == total_score:
                return True
        return False

    possible_bullseye_scores = set()

    # We use the simplified equation from Anna's score: 3*k1 + k2 + kb = 25,
    # where scores are s_i = 5 * k_i and k1 < k2 < k3 < kb.
    # k_i must be positive integers.

    # Iterate through all possible integer values for k1.
    # The loop can stop when 3*k1 is already greater than 25.
    for k1 in range(1, 9):
        # Iterate through all possible integer values for k2.
        # k2 must be greater than k1.
        for k2 in range(k1 + 1, 25):
            # Calculate kb based on k1 and k2.
            kb = 25 - (3 * k1) - k2

            # If the calculated kb is not positive, we can stop this inner loop.
            if kb <= 0:
                break

            # The condition k1 < k2 < k3 < kb implies there must be space for k3.
            # This means kb must be at least k2 + 2.
            if kb > k2 + 1:
                # We have a candidate set of scores (s1, s2, sb).
                s1 = 5 * k1
                s2 = 5 * k2
                sb = 5 * kb

                # Now, check if there exists a valid s3 between s2 and sb.
                for k3 in range(k2 + 1, kb):
                    s3 = 5 * k3
                    current_scores = [s1, s2, s3, sb]

                    # Check if this set of scores works for Bobby (total 230).
                    bobby_possible = is_sum_possible(current_scores, 230, 5)
                    
                    # If Bobby's score is possible, check Cliff's score (total 185).
                    if bobby_possible:
                        cliff_possible = is_sum_possible(current_scores, 185, 5)

                        # If both are possible, this sb is a valid bull's eye score.
                        if cliff_possible:
                            possible_bullseye_scores.add(sb)
                            # We only need one valid s3, so we can break and
                            # move to the next candidate for (s1, s2, sb).
                            break
    
    # Print the final answer.
    print(f"The exact number of possible values for the score of the bull's eye is: {len(possible_bullseye_scores)}")

# Execute the function to find the solution.
solve_archery_challenge()