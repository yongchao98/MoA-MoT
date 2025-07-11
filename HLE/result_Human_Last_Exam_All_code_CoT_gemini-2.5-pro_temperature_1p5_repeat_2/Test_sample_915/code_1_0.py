def is_score_possible(total_score, scores, num_arrows=5):
    """
    Checks if a total_score can be achieved by shooting num_arrows at a target with given scores.
    Iterates through all combinations of placing num_arrows in 4 scoring zones.
    """
    s_or, s_r1, s_r2, s_be = scores
    
    for n_be in range(num_arrows + 1):
        score_be = n_be * s_be
        if score_be > total_score:
            continue
            
        for n_r2 in range(num_arrows - n_be + 1):
            score_be_r2 = score_be + n_r2 * s_r2
            if score_be_r2 > total_score:
                continue

            for n_r1 in range(num_arrows - n_be - n_r2 + 1):
                n_or = num_arrows - n_be - n_r2 - n_r1
                
                current_score = (score_be_r2 + 
                                 n_r1 * s_r1 + 
                                 n_or * s_or)
                
                if current_score == total_score:
                    return True
    return False

def find_bulls_eye_values():
    """
    Finds all possible values for the bull's eye score based on the problem's constraints.
    """
    # Using a set to store solutions to only print unique s_be values and their first found valid equation
    found_solutions = set()
    
    # Iterate through possible values for 'a', which corresponds to s_or = 5*a
    # Constraint from Cliff's score: 5*s_or <= 185 => s_or <= 37 => 5a <= 37 => a <= 7
    # Constraint from combining Anna's equation and d>=10: 4a+2 <= 15 => 4a <= 13 => a <= 3
    for a in range(1, 4):
        s_or = 5 * a
        
        # Iterate through possible values for 'c', which corresponds to s_r2 = 5*c
        # Constraint from score ordering: c >= a + 2
        # Constraint from Anna's equation and d>=10: c <= 15 - 3a
        # Constraint from Anna's equation and d>c: c < (25 - 3a)/2
        c_min = a + 2
        c_max = min(15 - 3 * a, (25 - 3 * a - 1) // 2)
        
        for c in range(c_min, c_max + 1):
            s_r2 = 5 * c
            
            # From Anna's equation (3a + c + d = 25), calculate d
            d = 25 - 3 * a - c
            s_be = 5 * d
            
            # Constraint from Bobby's score: 5*s_be >= 230 => s_be >= 46 => d >= 10.
            # This is already incorporated into the loop range for 'c'.
            
            # If we've already found a valid scenario for this s_be, skip.
            if s_be in {sol[0] for sol in found_solutions}:
                continue

            # Iterate through possible values for 'b' (s_r1) where a < b < c
            for b in range(a + 1, c):
                s_r1 = 5 * b
                
                scores = (s_or, s_r1, s_r2, s_be)
                
                # A score set is valid if it allows for both Bobby's and Cliff's scores.
                bobby_possible = is_score_possible(230, scores)
                if bobby_possible:
                    cliff_possible = is_score_possible(185, scores)
                    if cliff_possible:
                        # Found a valid score set. Add the solution for printing and break the 'b' loop.
                        found_solutions.add((s_be, s_or, s_r2))
                        break

    # Print the results
    if not found_solutions:
        print("No possible values found for the bull's eye score.")
    else:
        # Sort solutions by bull's eye score for clean output
        sorted_solutions = sorted(list(found_solutions))
        for s_be, s_or, s_r2 in sorted_solutions:
            print(f"Possible bull's eye score: {s_be}. Validated with Anna's equation: 3 * {s_or} + {s_r2} + {s_be} = {3 * s_or + s_r2 + s_be}")
        
    print(f"\nThe exact number of possible values for the score of the bull's eye is: {len(found_solutions)}")
    
    # Final answer format
    print(f"\n<<<{len(found_solutions)}>>>")

if __name__ == '__main__':
    find_bulls_eye_values()