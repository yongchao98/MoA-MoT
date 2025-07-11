import sys

def check_score(total, s1, s2, s3, s4):
    """
    Checks if a total score can be achieved with 5 arrows hitting targets
    with scores {s1, s2, s3, s4}.
    """
    for n4 in range(6):  # Number of hits on the bull's eye
        for n3 in range(6 - n4):
            for n2 in range(6 - n4 - n3):
                n1 = 5 - (n4 + n3 + n2)
                # n1 is determined by the other counts
                
                current_total = n1 * s1 + n2 * s2 + n3 * s3 + n4 * s4
                if current_total == total:
                    return True
    return False

def solve_archery_puzzle():
    """
    Finds the number of possible values for the bull's eye score.
    """
    possible_s4_values = set()
    bobby_score = 230
    cliff_score = 185
    anna_total = 125

    # s4 must be at least ceil(230/5) = 46, so >= 50 as it's a multiple of 5.
    # s4 must be less than 125.
    for s4 in range(50, 125, 5):
        is_s4_plausible = False
        
        # From Anna's score: 3*s1 + s2 + s4 = 125 => 3*s1 + s2 = 125 - s4
        # We also know s1 >= 5.
        # Loop through possible s1 values.
        max_s1 = (anna_total - s4) // 3
        for s1 in range(5, max_s1 + 1, 5):
            s2 = anna_total - s4 - 3 * s1

            # A valid scenario requires s1 < s2 < s4 and a possible s3.
            # A slot for s3 exists if s4 - s2 >= 10.
            if s1 < s2 < s4 and (s4 - s2) >= 10:
                
                # Now, check if there EXISTS an s3 that makes Bobby's and Cliff's scores possible.
                for s3 in range(s2 + 5, s4, 5):
                    
                    bobby_possible = check_score(bobby_score, s1, s2, s3, s4)
                    if bobby_possible:
                        cliff_possible = check_score(cliff_score, s1, s2, s3, s4)
                        if cliff_possible:
                            # This s4 value is possible because we found a valid set of scores.
                            is_s4_plausible = True
                            
                            # Print the details for the first valid scenario found for this s4
                            if s4 not in possible_s4_values:
                                print(f"Found a valid scenario for bull's eye score s4={s4}:")
                                print(f"  Possible scores: s1={s1}, s2={s2}, s3={s3}, s4={s4}")
                                print(f"  Anna's equation: 3 * {s1} + {s2} + {s4} = {3*s1 + s2 + s4}")
                                sys.stdout.flush()
                            break  # Found a valid s3, no need to check others for this s1/s2 pair
            
            if is_s4_plausible:
                break # Found a valid s1/s2/s3, no need to check other s1 values for this s4
                
        if is_s4_plausible:
            possible_s4_values.add(s4)
    
    print("\n-------------------------------------------------")
    print(f"The possible values for the bull's eye score are: {sorted(list(possible_s4_values))}")
    print(f"The total number of possible values is: {len(possible_s4_values)}")

solve_archery_puzzle()