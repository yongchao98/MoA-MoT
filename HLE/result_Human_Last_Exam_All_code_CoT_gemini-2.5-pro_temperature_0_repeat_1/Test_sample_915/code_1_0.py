import math

def can_achieve_score(s1, s2, s3, s4, total):
    """
    Checks if a total score can be achieved with 5 arrows on a target
    with scores s1, s2, s3, s4.
    """
    num_arrows = 5
    # n4, n3, n2, n1 are the number of hits on the bull's eye, inner ring, etc.
    for n4 in range(num_arrows + 1):
        for n3 in range(num_arrows - n4 + 1):
            for n2 in range(num_arrows - n4 - n3 + 1):
                n1 = num_arrows - n4 - n3 - n2
                current_score = n1 * s1 + n2 * s2 + n3 * s3 + n4 * s4
                if current_score == total:
                    return True
    return False

def solve_archery_puzzle():
    """
    Finds the number of possible values for the bull's eye score.
    """
    # A set to store unique possible values for the bull's eye score (s4)
    possible_bullseye_scores = set()

    # Anna's equation: 3*s1 + s2 + s4 = 125
    # Scores are multiples of 5, and 0 < s1 < s2 < s3 < s4.

    # From 3*s1 + s2 + s4 = 125 and s1 < s2 < s4, we can deduce
    # 3*s1 + (s1+5) + (s1+10) <= 125 => 5*s1 <= 110 => s1 <= 22.
    # So, s1 can be 5, 10, 15, 20.
    for s1 in range(5, 25, 5):
        
        # From 3*s1 + s2 + s4 = 125 and s2 < s4, we get 2*s2 < 125 - 3*s1.
        s2_limit = math.ceil((125 - 3 * s1) / 2)
        for s2 in range(s1 + 5, s2_limit, 5):
            
            # Calculate s4 based on Anna's score
            s4 = 125 - 3 * s1 - s2

            # The scores must be strictly increasing (s2 < s4)
            if s4 <= s2:
                continue
            
            # Optimization: Bobby's score is 230. The max he can get is 5*s4.
            # So, 5*s4 >= 230 => s4 >= 46. Since s4 is a multiple of 5, s4 >= 50.
            if s4 < 50:
                continue

            # Iterate through possible values for s3 (s2 < s3 < s4)
            for s3 in range(s2 + 5, s4, 5):
                
                # Check if this score set is possible for Bobby (score 230)
                if not can_achieve_score(s1, s2, s3, s4, 230):
                    continue

                # Check if this score set is possible for Cliff (score 185)
                if not can_achieve_score(s1, s2, s3, s4, 185):
                    continue
                
                # If both checks pass, this is a valid scenario.
                # Add the bull's eye score to our set of possibilities.
                possible_bullseye_scores.add(s4)

    # The final answer is the number of unique values found for s4.
    final_count = len(possible_bullseye_scores)
    
    print("The possible values for the bull's eye score are:", sorted(list(possible_bullseye_scores)))
    print("The final equation is: Number of possible values = " + str(final_count))

solve_archery_puzzle()