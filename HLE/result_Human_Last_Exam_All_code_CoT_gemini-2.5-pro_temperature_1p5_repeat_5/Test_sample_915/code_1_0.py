import math

def solve():
    """
    Finds the number of possible values for the bull's eye score by systematically
    checking all combinations of scores that satisfy the problem's constraints.
    """

    def can_make_score(target, num_arrows, scores):
        """
        Checks if a target score can be achieved by summing num_arrows shots,
        with each shot having a score from the 'scores' tuple.
        """
        s1, s2, s3, sb = scores
        # Iterate through all combinations of 5 arrows for the 4 target areas
        for num_sb in range(num_arrows + 1):
            for num_s3 in range(num_arrows - num_sb + 1):
                for num_s2 in range(num_arrows - num_sb - num_s3 + 1):
                    num_s1 = num_arrows - num_sb - num_s3 - num_s2
                    # This check is implicitly handled by the loop ranges, but good practice.
                    if num_s1 >= 0:
                        calculated_score = (num_s1 * s1) + (num_s2 * s2) + (num_s3 * s3) + (num_sb * sb)
                        if calculated_score == target:
                            return True
        return False

    possible_bullseye_scores = set()

    # The scores s1, s2, s3, sb are positive multiples of 5, and s1 < s2 < s3 < sb.
    # Anna's score equation: 3*s1 + s2 + sb = 125.
    
    # From inequalities, 125 = 3*s1 + s2 + sb > 3*s1 + s1 + s1 = 5*s1, which means s1 < 25.
    # As s1 is a positive multiple of 5, s1 can be 5, 10, 15, or 20.
    for s1 in range(5, 25, 5):
        
        # We can also constrain sb for a given s1.
        # sb > s2 => sb > 125 - 3*s1 - sb => 2*sb > 125 - 3*s1
        # s2 > s1 => 125 - 3*s1 - sb > s1 => 125 - 4*s1 > sb
        sb_lower_bound = (125 - 3 * s1) / 2
        sb_upper_bound = 125 - 4 * s1
        
        # Iterate through possible bull's eye scores (sb), which must be multiples of 5.
        # Calculate the smallest multiple of 5 that is strictly greater than the lower bound.
        start_sb = math.ceil((sb_lower_bound + 1) / 5) * 5
        
        for sb in range(start_sb, sb_upper_bound, 5):
            
            # Calculate s2 based on Anna's equation.
            s2 = 125 - 3 * s1 - sb

            # A valid s3 must exist between s2 and sb. This requires sb and s2
            # to be separated by more than 5.
            if sb - s2 <= 5:
                continue

            # Iterate through possible s3 values.
            for s3 in range(s2 + 5, sb, 5):
                
                scores = (s1, s2, s3, sb)

                # Check if this set of scores works for Bobby (230) and Cliff (185).
                bobby_ok = can_make_score(230, 5, scores)
                if bobby_ok:
                    cliff_ok = can_make_score(185, 5, scores)
                    if cliff_ok:
                        # If a valid set of scores is found, the bull's eye score sb is possible.
                        possible_bullseye_scores.add(sb)
                        # We only need to find one valid s3 for a given sb.
                        # We can break and check the next potential sb.
                        break 
                        
    # The question asks for the number of possible values for the bull's eye.
    # The final equation is effectively the count of the valid possibilities found.
    print(len(possible_bullseye_scores))

solve()