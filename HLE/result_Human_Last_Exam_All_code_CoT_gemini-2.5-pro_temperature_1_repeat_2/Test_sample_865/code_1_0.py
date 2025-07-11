def solve_crease_pattern(patterns):
    """
    Analyzes a list of single-vertex crease patterns to find the angle 't'
    that makes them flat-foldable.

    A pattern is flat-foldable if it satisfies:
    1. Maekawa's Theorem: |Number of Mountains - Number of Valleys| = 2
    2. Kawasaki's Theorem: The sum of alternating angles is 180 degrees.
    """
    final_results = []
    
    for i, pattern_data in enumerate(patterns):
        pattern_num = i + 1
        print(f"--- Analyzing Pattern {pattern_num}: {pattern_data} ---")
        
        # 1. Parse the input pattern
        angles_str = [pattern_data[j] for j in range(0, len(pattern_data), 2)]
        creases = [pattern_data[j] for j in range(1, len(pattern_data), 2)]
        
        # 2. Check Maekawa's Theorem
        num_m = creases.count('M')
        num_v = creases.count('V')
        
        print(f"Checking Maekawa's Theorem...")
        print(f"Number of mountain (M) folds: {num_m}")
        print(f"Number of valley (V) folds: {num_v}")
        
        if abs(num_m - num_v) != 2:
            print(f"|M - V| = |{num_m} - {num_v}| = {abs(num_m - num_v)}. This must be 2.")
            print("Result for Pattern {}: none\n".format(pattern_num))
            final_results.append('none')
            continue
        
        print(f"|M - V| = |{num_m} - {num_v}| = 2. Maekawa's Theorem is satisfied.")
        
        # 3. Check Kawasaki's Theorem and solve for 't'
        print("Checking Kawasaki's Theorem...")
        
        angles = []
        t_is_present = 't' in angles_str
        t_index = -1
        if t_is_present:
            t_index = angles_str.index('t')
            
        for a in angles_str:
            angles.append(a if a == 't' else int(a))
            
        odd_angles_all = [angles[j] for j in range(0, len(angles), 2)]
        even_angles_all = [angles[j] for j in range(1, len(angles), 2)]
        
        t_is_in_odd = t_is_present and (t_index % 2 == 0)
        
        solution = 'none'

        if t_is_in_odd:
            # 't' is in the odd-indexed set, so check the even-indexed set first
            known_angles = even_angles_all
            angles_to_solve = odd_angles_all
            known_sum = sum(known_angles)
            equation_str = " + ".join(map(str, known_angles))
            print(f"Sum of even-indexed angles: {equation_str} = {known_sum}")
            
            if known_sum == 180:
                # Now solve for 't' in the other set
                angles_to_solve.remove('t')
                known_part_of_solve_sum = sum(angles_to_solve)
                t = 180 - known_part_of_solve_sum
                
                equation_str = " + ".join(map(str, angles_to_solve)) + " + t"
                print(f"Sum of odd-indexed angles must also be 180: {equation_str} = 180")
                print(f"Solving for t: {known_part_of_solve_sum} + t = 180  =>  t = {t}")
                
                if t > 0:
                    solution = t
                else:
                    print(f"Calculated angle t={t} is not positive.")
            else:
                print("Sum is not 180. Kawasaki's Theorem fails.")
        
        else: # 't' is in the even-indexed set or not present
            # Check the odd-indexed set first
            known_angles = odd_angles_all
            angles_to_solve = even_angles_all
            known_sum = sum(known_angles)
            equation_str = " + ".join(map(str, known_angles))
            print(f"Sum of odd-indexed angles: {equation_str} = {known_sum}")
            
            if known_sum == 180:
                if t_is_present: # 't' is in even set
                    angles_to_solve.remove('t')
                    known_part_of_solve_sum = sum(angles_to_solve)
                    t = 180 - known_part_of_solve_sum
                    
                    equation_str = " + ".join(map(str, angles_to_solve)) + " + t"
                    print(f"Sum of even-indexed angles must also be 180: {equation_str} = 180")
                    print(f"Solving for t: {known_part_of_solve_sum} + t = 180  =>  t = {t}")

                    if t > 0:
                        solution = t
                    else:
                        print(f"Calculated angle t={t} is not positive.")
                else: # 't' is not present, just check if valid
                    even_sum = sum(even_angles_all)
                    equation_str_even = " + ".join(map(str, even_angles_all))
                    print(f"Sum of even-indexed angles: {equation_str_even} = {even_sum}")
                    if even_sum != 180:
                        print("Sum is not 180. Kawasaki's Theorem fails.")
                    else: # Valid pattern, but no 't' to solve for.
                        print("Pattern is valid, but no 't' to solve for.")

            else:
                print("Sum is not 180. Kawasaki's Theorem fails.")
        
        print(f"Result for Pattern {pattern_num}: {solution}\n")
        final_results.append(solution)
        
    return final_results

if __name__ == '__main__':
    patterns = [
        [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
        [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
        [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
        [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
    ]
    
    results = solve_crease_pattern(patterns)
    
    # This line is for the final answer format requested by the user prompt
    # and would not typically be part of a general-purpose script.
    print("--- Summary ---")
    print(f"Final comma-separated list of values: [{','.join(map(str, results))}]")
