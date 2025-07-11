def is_achievable(target, num_arrows, scores):
    """
    Checks if a target score can be achieved by summing num_arrows scores
    from the provided list of possible scores.
    Iterates through all combinations of 5 shots.
    """
    s1, s2, s3, s4 = scores
    # Quick checks for optimization
    if target % 5 != 0:
        return False
    if not (num_arrows * s1 <= target <= num_arrows * s4):
        return False

    for n1 in range(num_arrows + 1):
        for n2 in range(num_arrows - n1 + 1):
            for n3 in range(num_arrows - n1 - n2 + 1):
                n4 = num_arrows - n1 - n2 - n3
                
                # Calculate the score for this combination of hits
                current_score = n1 * s1 + n2 * s2 + n3 * s3 + n4 * s4
                
                if current_score == target:
                    return True
    return False

def find_possible_bullseye_scores():
    """
    Calculates the number of possible values for the bull's eye score based on all constraints.
    """
    possible_s4_values = set()
    bobby_score = 230
    cliff_score = 185
    num_arrows = 5
    
    # 1. Generate candidate (k1, k2, k4) tuples based on combined constraints.
    # Anna: 3*k1 + k2 + k4 = 25
    # Order: k1 < k2 and k4 >= k2 + 2
    # Bobby/Cliff: k1 <= 7 and k4 >= 10
    candidate_k_tuples = []
    for k1 in range(1, 5):  # Anna's equation implies 5*k1 < 23, so k1 < 4.6
        # Determine the valid range for k2
        k2_max_from_order = (23 - 3 * k1) // 2  # from k4 >= k2 + 2
        k2_max_from_min_k4 = 15 - 3 * k1       # from k4 >= 10
        k2_max = min(k2_max_from_order, k2_max_from_min_k4)
        
        for k2 in range(k1 + 1, k2_max + 1):
            k4 = 25 - 3 * k1 - k2
            candidate_k_tuples.append((k1, k2, k4))
            
    # 2. For each candidate (k1, k2, k4), check if a valid k3 exists.
    for k1, k2, k4 in candidate_k_tuples:
        s4 = 5 * k4
        # Iterate through all possible k3 values to find one that works.
        for k3 in range(k2 + 1, k4):
            s1 = 5 * k1
            s2 = 5 * k2
            s3 = 5 * k3
            
            scores = (s1, s2, s3, s4)
            
            # Check if this set of scores works for both Bobby and Cliff
            bobby_ok = is_achievable(bobby_score, num_arrows, scores)
            if not bobby_ok:
                continue  # Try the next k3
            
            cliff_ok = is_achievable(cliff_score, num_arrows, scores)
            
            # If both are achievable, this s4 is a possible value.
            if cliff_ok:
                possible_s4_values.add(s4)
                # Found a working scenario, no need to check other k3 for this (k1,k2,k4)
                break 

    print(f"The possible values for the bull's eye score are: {sorted(list(possible_s4_values))}")
    print(f"The number of possible values is: {len(possible_s4_values)}")

# Run the analysis
find_possible_bullseye_scores()