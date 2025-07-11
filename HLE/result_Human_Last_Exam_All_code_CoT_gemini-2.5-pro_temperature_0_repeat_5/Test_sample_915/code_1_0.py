import math
from itertools import combinations_with_replacement

def is_sum_possible(target_score, num_arrows, score_values):
    """
    Checks if a target score can be achieved by summing a specific number of scores
    from a given list of score values. It uses combinations_with_replacement to
    check all possible ways the arrows could have landed.
    """
    for combo in combinations_with_replacement(score_values, num_arrows):
        if sum(combo) == target_score:
            return True
    return False

def solve_archery_puzzle():
    """
    Solves the archery puzzle to find the number of possible values for the bull's eye.

    The function iterates through all possible score sets that satisfy Anna's condition
    and then validates them against Bobby's and Cliff's scores.
    """
    # Let the scores be s1, s2, s3, s4. They are positive multiples of 5.
    # Let s_i = 5 * k_i, where k_i are positive integers with k1 < k2 < k3 < k4.
    # Anna's score: 3*s1 + s2 + s4 = 125  =>  3*k1 + k2 + k4 = 25.
    
    possible_s4_values = set()

    # From 3*k1 + k2 + k4 = 25 and k1 < k2 < k4, we can deduce that k1 must be in {1, 2, 3}.
    # If k1=4, s1=20. Then s2+s4=65. Smallest s2 is 25, so s4=40.
    # The scores would be {20, 25, s3, 40}. The maximum possible score for 5 arrows is 5*40=200.
    # This contradicts Bobby's score of 230. So, k1 cannot be 4.
    for k1 in range(1, 4):
        s1 = 5 * k1
        
        # From k2 < k4 and k2 + k4 = 25 - 3*k1, we get 2*k2 < 25 - 3*k1.
        k2_upper_bound = math.ceil((25 - 3 * k1) / 2)

        for k2 in range(k1 + 1, k2_upper_bound):
            s2 = 5 * k2
            k4 = 25 - 3 * k1 - k2
            s4 = 5 * k4

            # For a valid s3 to exist, we need k4 > k2 + 1.
            if k4 > k2 + 1:
                # We only need to find one valid s3 for this s4 to be possible.
                for k3 in range(k2 + 1, k4):
                    s3 = 5 * k3
                    scores = [s1, s2, s3, s4]

                    # Check if this set of scores works for Bobby (230) and Cliff (185).
                    bobby_possible = is_sum_possible(230, 5, scores)
                    if bobby_possible:
                        cliff_possible = is_sum_possible(185, 5, scores)
                        if cliff_possible:
                            # A valid score set is found. Add s4 to our set of possibilities.
                            possible_s4_values.add(s4)
                            # Break the inner loop since we've confirmed this s4 is possible.
                            break
    
    print(f"The possible values for the bull's eye score are: {sorted(list(possible_s4_values))}")
    print(f"The number of possible values for the bull's eye is: {len(possible_s4_values)}")

solve_archery_puzzle()
<<<12>>>