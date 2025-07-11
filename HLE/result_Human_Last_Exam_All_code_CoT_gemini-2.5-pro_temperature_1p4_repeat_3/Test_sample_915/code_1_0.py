from itertools import combinations_with_replacement

def is_score_possible(scores, num_arrows, target_sum):
    """
    Checks if a target sum can be achieved by summing a number of scores
    chosen with replacement from the provided list of scores.
    """
    for combo in combinations_with_replacement(scores, num_arrows):
        if sum(combo) == target_sum:
            return True
    return False

def solve_archery_puzzle():
    """
    Finds the number of possible values for the bull's eye score based on
    the given information.
    """
    possible_bullseye_solutions = {}
    bobby_target = 230 // 5
    cliff_target = 185 // 5
    anna_target = 125 // 5
    
    # Iterate through possible values for x1 (outer ring score / 5)
    # x1 must be less than 4, otherwise Anna's score becomes too high.
    for x1 in range(1, 4):
        # Iterate through possible values for x2 (next ring score / 5)
        # We know x2 > x1 and from Anna's eq: x2 < x4. This implies 2*x2 < 25 - 3*x1
        for x2 in range(x1 + 1, (anna_target - 3 * x1 + 1) // 2):
            # Calculate x4 from Anna's score
            x4 = anna_target - 3 * x1 - x2

            # The score for the bull's eye must be strictly greater than for the inner ring (x3),
            # which must be strictly greater than x2. So there must be space for x3.
            if x4 <= x2 + 1:
                continue

            # Iterate through possible values for x3
            for x3 in range(x2 + 1, x4):
                scores = (x1, x2, x3, x4)
                
                # Check if this set of scores can explain Bobby's and Cliff's totals
                bobby_possible = is_score_possible(scores, 5, bobby_target)
                if not bobby_possible:
                    continue
                
                cliff_possible = is_score_possible(scores, 5, cliff_target)
                if not cliff_possible:
                    continue

                # If we reach here, we found a valid set of scores.
                # The bull's eye score is s4 = x4 * 5.
                s4 = x4 * 5
                
                # Store the solution and the score set that proved it.
                if s4 not in possible_bullseye_solutions:
                    s1 = x1 * 5
                    s2 = x2 * 5
                    s3 = x3 * 5
                    possible_bullseye_solutions[s4] = {
                        "s1": s1, "s2": s2, "s3": s3, "s4": s4
                    }
                
                # Once we find one valid x3 for this (x1, x2, x4) triplet,
                # this bull's eye value is confirmed as possible, so we can stop checking other x3s.
                break
    
    # Print the final results
    num_solutions = len(possible_bullseye_solutions)
    print(f"Found {num_solutions} possible value(s) for the bull's eye score.\n")
    
    sorted_solutions = sorted(possible_bullseye_solutions.keys())

    for s4 in sorted_solutions:
        solution_set = possible_bullseye_solutions[s4]
        s1, s2 = solution_set['s1'], solution_set['s2']
        
        print(f"Possible Bull's Eye Score: {s4}")
        print(f"  - An example of a valid full score set is: {tuple(solution_set.values())}")
        print(f"  - Verifying with Anna's score: 3 * {s1} + 1 * {s2} + 1 * {s4} = {3*s1 + s2 + s4} (which is 125 as required).")
        print("-" * 20)

    print(f"\nBased on this analysis, the exact number of possible values for the bull's eye is {num_solutions}.")
    print(f"<<<{num_solutions}>>>")

# Run the solver
solve_archery_puzzle()