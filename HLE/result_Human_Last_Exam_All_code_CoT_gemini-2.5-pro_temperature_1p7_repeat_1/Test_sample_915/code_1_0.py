import itertools

def is_achievable(target_sum, scores, num_shots=5):
    """
    Checks if a target sum can be achieved by summing num_shots scores
    from the given list of scores.
    """
    for combo in itertools.combinations_with_replacement(scores, num_shots):
        if sum(combo) == target_sum:
            return True
    return False

def find_possible_bullseye_scores():
    """
    Finds the number of possible values for the bull's eye score based on the
    problem's constraints.
    """
    possible_bullseye_scores = set()
    anna_target_scaled = 125 // 5
    bobby_target_scaled = 230 // 5
    cliff_target_scaled = 185 // 5

    # From 3*x1 + x2 + x_bull = 25 and x1 < x2 < x3 < x_bull, we can
    # derive 5*x1 + 4 <= 25, which means x1 <= 4.2. So x1 is in {1, 2, 3, 4}.
    for x1 in range(1, 5):
        # We need x2 > x1
        # Loop for x2, the upper bound is determined by the equation and x_bull > x2
        for x2 in range(x1 + 1, anna_target_scaled):
            # From Anna's equation: 3*x1 + x2 + x_bull = 25
            x_bull = anna_target_scaled - 3 * x1 - x2

            # The scores must be strictly increasing. If x_bull is not greater
            # than x2, then increasing x2 further will only make x_bull smaller.
            # So, we can break from this inner loop.
            if x_bull <= x2:
                break

            # To fit x3 between x2 and x_bull (x2 < x3 < x_bull),
            # the gap between them must be greater than 1.
            if x_bull - x2 <= 1:
                continue

            # Iterate through all possible integer values for x3
            for x3 in range(x2 + 1, x_bull):
                scores_scaled = (x1, x2, x3, x_bull)

                # Check if this score set is valid for Bobby
                if not is_achievable(bobby_target_scaled, scores_scaled):
                    continue

                # Check if this score set is valid for Cliff
                if not is_achievable(cliff_target_scaled, scores_scaled):
                    continue

                # If the score set is valid for everyone, add the actual
                # bull's eye score to our set of possibilities.
                possible_bullseye_scores.add(x_bull * 5)
    
    return possible_bullseye_scores

# Run the solver
final_scores = find_possible_bullseye_scores()
final_count = len(final_scores)

# Print the results as per the specified format
print(f"The possible values for the bull's eye score are: {sorted(list(final_scores))}")
print("The final equation for the total count is:")
# Output each number in the 'final equation' as requested.
print(f"Number of possible values = {final_count}")