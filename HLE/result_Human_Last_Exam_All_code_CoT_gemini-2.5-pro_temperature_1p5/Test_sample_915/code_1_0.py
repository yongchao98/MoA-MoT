import itertools

def is_achievable(scores, target_score, num_arrows=5):
    """
    Checks if a target_score can be achieved by summing num_arrows values from the scores list.
    This is done by checking all combinations with replacement.
    """
    for hits in itertools.combinations_with_replacement(scores, num_arrows):
        if sum(hits) == target_score:
            return True
    return False

def solve_archery_puzzle():
    """
    Solves the archery puzzle to find the number of possible bull's eye scores.
    """
    # From Anna: 3*s1 + s2 + s4 = 125. Let s_i = 5*k_i.
    # -> 3*k1 + k2 + k4 = 25.
    # Conditions: 0 < k1 < k2 < s3/5 < k4.
    
    anna_triplets = []
    # k1 can be at most 4, otherwise 3*k1 + (k1+1) + (k1+2) > 25.
    for k1 in range(1, 5):
        # k2 must be greater than k1. The upper bound is loose.
        for k2 in range(k1 + 1, 25):
            k4 = 25 - 3 * k1 - k2
            # Check conditions: 
            # 1. k4 must be an integer > k2.
            # 2. There must be space for k3, so k4 > k2 + 1.
            if k4 > k2 + 1:
                anna_triplets.append((k1, k2, k4))

    possible_s4_values = set()

    # For each valid triplet from Anna's score, check Bobby and Cliff.
    for k1, k2, k4 in anna_triplets:
        s1 = 5 * k1
        s2 = 5 * k2
        s4 = 5 * k4

        # Iterate through possible s3 values.
        # k3 must be between k2 and k4.
        for k3 in range(k2 + 1, k4):
            s3 = 5 * k3
            scores = (s1, s2, s3, s4)
            
            # Check if this score set allows Bobby and Cliff to get their scores.
            bobby_possible = is_achievable(scores, 230)
            cliff_possible = is_achievable(scores, 185)

            if bobby_possible and cliff_possible:
                # If a valid set of scores is found, this s4 is possible.
                # We can add it to our set and stop checking other s3 values
                # for this (s1, s2, s4) combination.
                possible_s4_values.add(s4)
                break
    
    # The final answer is the number of unique possible values for the bull's eye.
    final_count = len(possible_s4_values)
    print(f"The exact number of possible values for the score of the bull's eye is {final_count}.")

solve_archery_puzzle()