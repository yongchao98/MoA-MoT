import sys

def is_score_possible(s1, s2, s3, s4, total_arrows, total_score):
    """
    Checks if a total score can be achieved by iterating through all 
    combinations of 5 shots across the four scoring regions.
    n1, n2, n3, n4 are the number of hits on each respective ring.
    """
    for n4 in range(total_arrows + 1):
        for n3 in range(total_arrows - n4 + 1):
            for n2 in range(total_arrows - n4 - n3 + 1):
                n1 = total_arrows - n4 - n3 - n2
                
                current_score = n1 * s1 + n2 * s2 + n3 * s3 + n4 * s4
                if current_score == total_score:
                    return True
    return False

def find_bullseye_possibilities():
    """
    Finds all possible values for the bull's eye score (s4) based on the
    information given for Anna, Bobby, and Cliff.
    """
    # A set to store unique possible values for s4.
    possible_bullseye_scores = set()

    # Scores are s1 < s2 < s3 < s4, all multiples of 5.
    # Anna: 3*s1 + s2 + s4 = 125.
    # From constraints, we can deduce s1 can only be 5, 10, 15, or 20.
    for s1 in range(5, 25, 5):
        
        # From 3*s1 + s2 + s4 = 125 and s4 >= s2 + 10, we can determine the
        # valid range for s2.
        # 2*s2 <= 115 - 3*s1
        s2_upper_bound = (115 - 3 * s1) // 2

        for s2 in range(s1 + 5, s2_upper_bound + 1, 5):
            # s4 is determined by Anna's score.
            s4 = 125 - 3 * s1 - s2

            # Now, check if there exists a valid s3 between s2 and s4.
            found_valid_s3 = False
            for s3 in range(s2 + 5, s4, 5):
                # We have a full candidate score set (s1, s2, s3, s4).
                # Check if it works for both Bobby and Cliff.
                bobby_ok = is_score_possible(s1, s2, s3, s4, 5, 230)
                if not bobby_ok:
                    continue

                cliff_ok = is_score_possible(s1, s2, s3, s4, 5, 185)
                
                # If a valid score system is found for this (s1,s2,s4) triplet,
                # we can stop checking other s3's and mark this s4 as possible.
                if cliff_ok:
                    found_valid_s3 = True
                    break
            
            if found_valid_s3:
                possible_bullseye_scores.add(s4)

    # Output the results
    sorted_scores = sorted(list(possible_bullseye_scores))
    print("The possible values for the bull's eye score are:")
    print(sorted_scores)
    print(f"\nThe number of possible values for the score of the bull's eye is: {len(sorted_scores)}")

if __name__ == '__main__':
    find_bullseye_possibilities()
