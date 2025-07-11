def solve_flat_foldable_angle(patterns):
    """
    Analyzes a list of single-vertex crease patterns to find the angle 't'
    that makes them flat-foldable.

    Args:
        patterns: A list of lists, where each inner list represents a crease pattern.
    """
    final_results = []
    print("Analyzing crease patterns based on Maekawa's and Kawasaki's theorems:")
    
    for i, pattern_list in enumerate(patterns):
        # Parse angles and folds
        angles = [item for j, item in enumerate(pattern_list) if j % 2 == 0]
        folds = [item for j, item in enumerate(pattern_list) if j % 2 != 0]

        # Check for even number of angles
        if len(angles) % 2 != 0:
            final_results.append('none')
            continue

        # 1. Check Maekawa's Theorem: |N_M - N_V| = 2
        num_m = folds.count('M')
        num_v = folds.count('V')
        if abs(num_m - num_v) != 2:
            final_results.append('none')
            continue

        # 2. Check Kawasaki's Theorem
        try:
            t_index = angles.index('t')
        except ValueError:
            # 't' is not in the pattern, and Maekawa's theorem failed for the given case.
            final_results.append('none')
            continue

        # Separate known angles into odd and even positions
        sum_odd_known = sum(angle for j, angle in enumerate(angles) if j % 2 == 0 and angle != 't')
        sum_even_known = sum(angle for j, angle in enumerate(angles) if j % 2 != 0 and angle != 't')

        solution = None
        # If t is at an odd-indexed position (1st, 3rd, etc.)
        if t_index % 2 == 0:
            # The sum of even-indexed angles must be 180
            if sum_even_known == 180:
                t = 180 - sum_odd_known
                if t > 0:
                    solution = t
        # If t is at an even-indexed position (2nd, 4th, etc.)
        else:
            # The sum of odd-indexed angles must be 180
            if sum_odd_known == 180:
                t = 180 - sum_even_known
                if t > 0:
                    solution = t
        
        if solution is not None:
            # Format the equation string for the successful case
            if t_index % 2 == 0: # t is in odd group
                equation_str = f"{' + '.join(map(str, [a for a in angles if angles.index(a) % 2 == 0 and a != 't']))} + t = 180"
            else: # t is in even group
                equation_str = f"{' + '.join(map(str, [a for a in angles if angles.index(a) % 2 != 0 and a != 't']))} + t = 180"
            
            solved_t_str = int(solution) if solution == int(solution) else solution
            print(f"For pattern {i+1}, solving the equation {equation_str} gives t = {solved_t_str}.")
            final_results.append(solved_t_str)
        else:
            final_results.append('none')

    # Print the final comma-separated list
    print("\nFinal comma-separated list of values for t:")
    print(f"[{','.join(map(str, final_results))}]")


# The four crease patterns from the problem
crease_patterns = [
    [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
    [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
    [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
    [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
]

solve_flat_foldable_angle(crease_patterns)