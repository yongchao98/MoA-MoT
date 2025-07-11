import sys

def solve_crease_pattern(pattern_name, pattern_data):
    """
    Analyzes a single-vertex crease pattern for flat-foldability and solves for 't'.
    Prints the analysis steps.
    """
    # Use a file-like object to capture print statements for the final equation
    from io import StringIO
    original_stdout = sys.stdout
    sys.stdout = captured_output = StringIO()

    try:
        angles = []
        folds = []
        t_present = False

        # 1. Parse the input into angles and folds
        for i, item in enumerate(pattern_data):
            if i % 2 == 0:  # Angle
                if item == 't':
                    angles.append('t')
                    t_present = True
                else:
                    angles.append(item)
            else:  # Fold type
                folds.append(item)

        # 2. Check Maekawa's Theorem: |M - V| = 2
        num_m = folds.count('M')
        num_v = folds.count('V')
        
        print(f"Analysis for pattern {pattern_name}: {pattern_data}")
        print(f"Step 1: Check Maekawa's Theorem (|M - V| = 2)")
        print(f"  - Mountain folds (M): {num_m}")
        print(f"  - Valley folds (V): {num_v}")
        if abs(num_m - num_v) != 2:
            print(f"  - Result: |{num_m} - {num_v}| = {abs(num_m-num_v)} != 2. The pattern is not flat-foldable.")
            return 'none'
        print(f"  - Result: |{num_m} - {num_v}| = 2. Maekawa's Theorem is satisfied.")

        # 3. Check Kawasaki's Theorem: Sum of alternating angles = 180
        print(f"\nStep 2: Check Kawasaki's Theorem (Alternating angle sums must be 180)")
        odd_angles = [angles[i] for i in range(0, len(angles), 2)]
        even_angles = [angles[i] for i in range(1, len(angles), 2)]
        
        # Case where there is no 't' to solve for
        if not t_present:
            sum_odd = sum(odd_angles)
            sum_even = sum(even_angles)
            print(f"  - Sum of odd-indexed angles: {' + '.join(map(str, odd_angles))} = {sum_odd}")
            print(f"  - Sum of even-indexed angles: {' + '.join(map(str, even_angles))} = {sum_even}")
            if sum_odd != 180 or sum_even != 180:
                print("  - Result: The sums are not both 180. The pattern is not flat-foldable.")
                return 'none'
            else:
                # This case isn't in the problem set, but would be valid if it occurred
                return 'valid but no t'

        # Case where 't' needs to be solved
        if 't' in odd_angles:
            # 't' is odd-indexed, so check the sum of even-indexed angles first
            sum_even = sum(even_angles)
            print(f"  - Checking the constant sum (even-indexed angles): {' + '.join(map(str, even_angles))} = {sum_even}")
            if sum_even != 180:
                print(f"  - Result: Sum is {sum_even}, not 180. The pattern is not flat-foldable.")
                return 'none'
            
            # Now solve for 't' using the odd-indexed angles
            known_odd = [a for a in odd_angles if a != 't']
            sum_known_odd = sum(known_odd)
            t_value = 180 - sum_known_odd
            print(f"  - Solving for 't' using the other sum (odd-indexed angles):")
            print(f"    {' + '.join(map(str, known_odd))} + t = 180")
            print(f"    {sum_known_odd} + t = 180")
            print(f"    t = 180 - {sum_known_odd}")

        else: # 't' is in even_angles
            # 't' is even-indexed, so check the sum of odd-indexed angles first
            sum_odd = sum(odd_angles)
            print(f"  - Checking the constant sum (odd-indexed angles): {' + '.join(map(str, odd_angles))} = {sum_odd}")
            if sum_odd != 180:
                print(f"  - Result: Sum is {sum_odd}, not 180. The pattern is not flat-foldable.")
                return 'none'
                
            # Now solve for 't' using the even-indexed angles
            known_even = [a for a in even_angles if a != 't']
            sum_known_even = sum(known_even)
            t_value = 180 - sum_known_even
            print(f"  - Solving for 't' using the other sum (even-indexed angles):")
            print(f"    {' + '.join(map(str, known_even))} + t = 180")
            print(f"    {sum_known_even} + t = 180")
            print(f"    t = 180 - {sum_known_even}")
        
        if t_value <= 0:
            print(f"  - Result: Calculated angle t = {t_value} must be positive. The pattern is not flat-foldable.")
            return 'none'

        print(f"  - Result: The required angle is t = {t_value}.")
        return t_value

    finally:
        # Restore stdout and print the captured output
        sys.stdout = original_stdout
        print(captured_output.getvalue())


def solve_all_patterns():
    patterns = {
        "1": [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
        "2": [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
        "3": [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
        "4": [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
    }

    results = []
    for name, data in patterns.items():
        result = solve_crease_pattern(name, data)
        results.append(result)
        print("-" * 30)

    # Print the final comma-separated list
    final_result_string = ",".join(map(str, results))
    print(f"\nFinal comma-separated list of values: [{final_result_string}]")
    
    # Required final answer format
    print(f"<<<[{final_result_string}]>>>")

if __name__ == '__main__':
    solve_all_patterns()