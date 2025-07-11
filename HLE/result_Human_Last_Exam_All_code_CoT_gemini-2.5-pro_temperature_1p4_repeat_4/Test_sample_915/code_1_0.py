import math

def check_total_score(target_k, k_values):
    """
    Checks if a target score (in k-units) can be formed by summing 5 scores from k_values.
    k_values = [k_r3, k_r2, k_r1, k_be]
    """
    num_shots = 5
    k_r3, k_r2, k_r1, k_be = k_values
    
    # n_be, n_r1, n_r2, n_r3 are the number of hits in each region
    for n_be in range(num_shots + 1):
        for n_r1 in range(num_shots - n_be + 1):
            for n_r2 in range(num_shots - n_be - n_r1 + 1):
                n_r3 = num_shots - n_be - n_r1 - n_r2
                
                current_sum_k = n_be * k_be + n_r1 * k_r1 + n_r2 * k_r2 + n_r3 * k_r3
                if current_sum_k == target_k:
                    return True
    return False

def solve():
    """
    Finds the number of possible values for the bull's eye score.
    """
    possible_solutions = {}
    
    TARGET_BOBBY_K = 230 // 5
    TARGET_CLIFF_K = 185 // 5
    
    # From problem analysis: 10 <= k_be < 20 and 1 <= k_r3 <= 4
    for k_be in range(10, 20):
        # We store the result for a given k_be once found.
        if k_be in possible_solutions:
            continue
            
        # A tighter bound derived from all inequalities is k_r3 <= 4.
        for k_r3 in range(1, 5):
            # From Anna's equation: k_r2 = 25 - k_be - 3*k_r3
            # Check if this combination is valid based on score ordering.
            # Constraints: (26 - 2*k_be)/3 < k_r3 < (25 - k_be)/4
            lower_bound_kr3 = (26.0 - 2.0 * k_be) / 3.0
            upper_bound_kr3 = (25.0 - k_be) / 4.0
            
            if k_r3 > lower_bound_kr3 and k_r3 < upper_bound_kr3:
                k_r2 = 25 - k_be - 3 * k_r3
                
                # Search for a valid k_r1 to complete the score set
                for k_r1 in range(k_r2 + 1, k_be):
                    k_values = [k_r3, k_r2, k_r1, k_be]
                    
                    # Check if this set of scores works for Cliff and Bobby
                    cliff_ok = check_total_score(TARGET_CLIFF_K, k_values)
                    if cliff_ok:
                        bobby_ok = check_total_score(TARGET_BOBBY_K, k_values)
                        if bobby_ok:
                            # Found a valid score set. Store it and break loops for this k_be.
                            possible_solutions[k_be] = (k_r3, k_r2, k_r1)
                            break
                
                if k_be in possible_solutions:
                    break
    
    print("Based on the problem information, we search for sets of scores {s_r3, s_r2, s_r1, s_be} that satisfy all conditions.")
    print("Here are the possible values for the bull's eye score, with an example of a valid score set for each:\n")

    sorted_k_be = sorted(possible_solutions.keys())

    for k_be in sorted_k_be:
        k_r3, k_r2, k_r1 = possible_solutions[k_be]
        s_be = k_be * 5
        s_r1 = k_r1 * 5
        s_r2 = k_r2 * 5
        s_r3 = k_r3 * 5
        
        print(f"Possible Bull's Eye Score: {s_be}")
        print(f"  - Example valid score set (s_r3, s_r2, s_r1, s_be): ({s_r3}, {s_r2}, {s_r1}, {s_be})")
        
        # Outputting the numbers in Anna's equation
        anna_sum = 3 * s_r3 + s_r2 + s_be
        print(f"  - Anna's equation check: 3 * {s_r3} + {s_r2} + {s_be} = {anna_sum}")
        print("-" * 30)

    num_possible_values = len(possible_solutions)
    print(f"\nThe exact number of possible values for the score of the bull's eye is: {num_possible_values}")

    return num_possible_values

# Execute the solver function
final_answer = solve()
print(f"<<<{final_answer}>>>")