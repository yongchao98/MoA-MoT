import itertools

def check_score_possible(target_sum, scores, num_arrows=5):
    """
    Checks if a target_sum can be achieved with a given number of arrows
    on a set of scores.
    
    Args:
        target_sum (int): The total score to check for.
        scores (tuple): A tuple of the four possible scores (k1, k2, k3, k4).
        num_arrows (int): The number of arrows shot.

    Returns:
        bool: True if the score is possible, False otherwise.
    """
    # Generate all possible combinations of scores for the given number of arrows
    for combo in itertools.combinations_with_replacement(scores, num_arrows):
        if sum(combo) == target_sum:
            return True
    return False

def find_possible_bullseye_scores():
    """
    Finds all possible values for the bull's eye score based on the
    problem's constraints.
    """
    # Use a set to store unique possible values for the bull's eye score (k4)
    possible_k4_values = set()
    
    # Simplified scores for Bobby and Cliff
    bobby_target = 230 // 5  # 46
    cliff_target = 185 // 5  # 37

    # Iterate through possible values for k1.
    # From 3*k1 + k2 + k4 = 25 and k1 < k2 < k4, we can deduce 5*k1 < 25, so k1 < 5.
    for k1 in range(1, 5):
        # Iterate through possible values for k4 based on constraints from Anna's equation.
        # k4 > k2 = 25 - 3*k1 - k4  =>  2*k4 > 25 - 3*k1  => k4 > (25 - 3*k1) / 2
        # k2 > k1 => 25 - 3*k1 - k4 > k1 => 25 > 4*k1 + k4 => k4 < 25 - 4*k1
        k4_lower_bound = int((25 - 3 * k1) / 2) + 1
        k4_upper_bound = 25 - 4 * k1
        
        for k4 in range(k4_lower_bound, k4_upper_bound):
            # Calculate k2 from Anna's equation
            k2 = 25 - 3 * k1 - k4
            
            # The loop ranges for k1 and k4 already ensure k1 < k2 < k4.
            # Now, iterate through possible values for k3.
            for k3 in range(k2 + 1, k4):
                scores = (k1, k2, k3, k4)
                
                # Check if this set of scores is valid for both Bobby and Cliff
                if check_score_possible(cliff_target, scores):
                    if check_score_possible(bobby_target, scores):
                        possible_k4_values.add(k4)

    # Print the results
    print("Found the following possible values for the bull's eye score:")
    # Sort the results for clear presentation
    sorted_scores = sorted([val * 5 for val in possible_k4_values])
    for score in sorted_scores:
        print(f"- {score}")
        
    print("\n---")
    print(f"The final equation is the count of these possibilities.")
    print(f"Number of possible values = {len(sorted_scores)}")

# Run the solver
find_possible_bullseye_scores()
<<<9>>>