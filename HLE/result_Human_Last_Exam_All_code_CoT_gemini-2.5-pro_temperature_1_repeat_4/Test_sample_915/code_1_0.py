import collections

def solve_archery_puzzle():
    """
    Solves the archery puzzle by iterating through all possible score combinations
    that satisfy the given constraints.
    """
    
    # This dictionary will store a valid (s1, s2) combination for each unique s_bull found.
    # We use this to print the example equations as requested.
    valid_combinations = collections.OrderedDict()

    # Derived constraints from the problem description
    ANNA_TOTAL_SCORE = 125
    # Cliff's score implies 5*s1 <= 185 -> s1 <= 37. As a multiple of 5, s1 <= 35.
    S1_MAX = 35
    # Bobby's score implies 5*s_bull >= 230 -> s_bull >= 46. As a multiple of 5, s_bull >= 50.
    S_BULL_MIN = 50

    print("Finding valid score combinations based on the rules...\n")
    
    # Iterate through possible values for s1 (outer ring score)
    for s1 in range(5, S1_MAX + 1, 5):
        
        # The condition s2 < s3 < s_bull means s_bull >= s2 + 10.
        # From Anna's equation: s_bull = 125 - 3*s1 - s2
        # So: 125 - 3*s1 - s2 >= s2 + 10  =>  115 - 3*s1 >= 2*s2
        s2_upper_bound = (115 - 3 * s1) // 2
        
        # Iterate through possible values for s2 (next ring score)
        for s2 in range(s1 + 5, s2_upper_bound + 1, 5):
            
            # Calculate s_bull based on Anna's total score
            s_bull = ANNA_TOTAL_SCORE - 3 * s1 - s2
            
            # Check if this combination meets all constraints
            if s_bull >= S_BULL_MIN and s_bull > s2:
                # If this s_bull value is new, store its combination
                if s_bull not in valid_combinations:
                    valid_combinations[s_bull] = (s1, s2)

    # Sort the found bull's eye scores for a clean presentation
    sorted_bull_scores = sorted(valid_combinations.keys())

    for s_bull in sorted_bull_scores:
        s1, s2 = valid_combinations[s_bull]
        print(f"Possible Bull's Eye Score: {s_bull}")
        print(f"This is possible if, for example, s1={s1} and s2={s2}.")
        # Outputting the numbers in Anna's equation
        print(f"Anna's score check: 3 * {s1} + {s2} + {s_bull} = {3*s1 + s2 + s_bull}\n")
    
    print("-" * 40)
    print(f"The set of all possible values for the bull's eye is:")
    print(set(sorted_bull_scores))
    print("\n" + "=" * 40)
    print(f"The exact number of possible values for the score of the bull's eye is: {len(valid_combinations)}")
    print("=" * 40)

solve_archery_puzzle()