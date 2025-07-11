def solve_crease_pattern(pattern_data, case_number):
    """
    Solves for the angle 't' that makes a single-vertex crease pattern flat-foldable.

    Args:
        pattern_data (list): A list representing the crease pattern, e.g., [100, 'M', 62, 'V', ...].
        case_number (int): The number of the case for printing.

    Returns:
        A string or integer representing the value of 't', or 'none' if no solution exists.
    """
    print(f"--- Case {case_number}: {pattern_data} ---")
    
    angles = []
    folds = []
    t_is_present = 't' in pattern_data

    for i in range(0, len(pattern_data), 2):
        angles.append(pattern_data[i])
        folds.append(pattern_data[i+1])

    # 1. Check Maekawa's Theorem: |#M - #V| = 2
    num_m = folds.count('M')
    num_v = folds.count('V')
    if abs(num_m - num_v) != 2:
        print(f"Result: No solution. The pattern fails Maekawa's Theorem.")
        print(f"Reason: The number of mountain folds is {num_m} and valley folds is {num_v}.")
        print(f"The difference |{num_m} - {num_v}| = {abs(num_m - num_v)}, which is not 2.\n")
        return 'none'

    # 2. Check Kawasaki's Theorem: Alternating angle sums must be 180 degrees.
    # This part is only reached if Maekawa's Theorem is satisfied.
    
    # If t is not in the pattern, we can only check if it's valid as-is.
    if not t_is_present:
         print("Result: No solution. The pattern is invalid as given and has no variable 't' to solve.")
         return 'none'

    alt_angles1 = []
    alt_angles2 = []
    
    # Separate angles into two alternating groups
    for i, angle in enumerate(angles):
        if i % 2 == 0:
            alt_angles1.append(angle)
        else:
            alt_angles2.append(angle)
            
    # Determine which group contains 't' and solve
    if 't' in alt_angles1:
        sum_known_alt2 = sum(alt_angles2)
        if sum_known_alt2 != 180:
            print(f"Result: No solution. The pattern fails Kawasaki's Theorem.")
            print(f"Reason: The sum of one set of alternating angles is {sum_known_alt2}, which is not 180.")
            return 'none'
        
        known_angles_in_t_group = [a for a in alt_angles1 if a != 't']
        sum_known_in_t_group = sum(known_angles_in_t_group)
        t = 180 - sum_known_in_t_group
        
    else: # 't' is in alt_angles2
        sum_known_alt1 = sum(alt_angles1)
        if sum_known_alt1 != 180:
            print(f"Result: No solution. The pattern fails Kawasaki's Theorem.")
            print(f"Reason: The sum of one set of alternating angles ({' + '.join(map(str, alt_angles1))}) is {sum_known_alt1}, which is not 180.\n")
            return 'none'

        known_angles_in_t_group = [a for a in alt_angles2 if a != 't']
        sum_known_in_t_group = sum(known_angles_in_t_group)
        t = 180 - sum_known_in_t_group

    if t > 0:
        t_val = int(t) if t.is_integer() else t
        print(f"Result: Solution found. t = {t_val}")
        print("Reason: Maekawa's theorem is satisfied. For Kawasaki's theorem, alternating angles must sum to 180.")
        
        # Display the equation that was solved for t
        equation_parts = ' + '.join(map(str, known_angles_in_t_group))
        final_equation_parts = ' + '.join(map(str, known_angles_in_t_group + [t_val]))
        print(f"The governing equation is: {equation_parts} + t = 180.")
        print(f"The final equation with the solution is: {final_equation_parts} = 180.\n")
        return t_val
    else:
        print(f"Result: No solution. Solving for 't' yields a non-positive angle ({t}).\n")
        return 'none'

if __name__ == '__main__':
    # The four crease patterns from the problem
    patterns = [
        [100,'M',62,'V',22,'M','t','V',33,'M',90,'V'],
        [90,'M',120,'M',60,'M',90,'M'],
        [60,'V',60,'M',120,'M','t','M'],
        [77,'M',15,'M',50,'V',33,'M','t','V',130,'M']
    ]

    solutions = []
    for i, p in enumerate(patterns):
        solution = solve_crease_pattern(p, i + 1)
        solutions.append(solution)
        
    print("-----------------------------------------")
    print(f"Final comma-separated list of values: {solutions}")