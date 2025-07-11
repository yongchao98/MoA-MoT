import math

def solve_crease_pattern(pattern_data, pattern_id):
    """
    Analyzes a single-vertex crease pattern to find the angle 't' that makes it flat-foldable.
    Prints the analysis steps.
    Returns the value of 't' or 'none'.
    """
    print(f"Analyzing pattern {pattern_id}: {pattern_data}")
    
    angles_str = []
    creases = []
    has_t = False
    
    # 1. Parse input pattern
    for i, item in enumerate(pattern_data):
        if i % 2 == 0:
            angles_str.append(item)
            if item == 't':
                has_t = True
        else:
            creases.append(item)

    # 2. Check Maekawa's Theorem
    num_m = creases.count('M')
    num_v = creases.count('V')
    print(f"Step 1: Maekawa's Theorem Check")
    print(f"  Mountain folds (M) = {num_m}, Valley folds (V) = {num_v}")
    if abs(num_m - num_v) != 2:
        print(f"  Check FAILED: |#M - #V| = |{num_m} - {num_v}| = {abs(num_m - num_v)}, which is not 2.")
        print(f"Result for pattern {pattern_id} is 'none'.\n")
        return "none"
    print(f"  Check PASSED: |#M - #V| = |{num_m} - {num_v}| = 2.")

    # Handle case with no 't' variable
    if not has_t:
        print("Step 2: Angle Check")
        print("  Pattern has no variable 't' to solve for.")
        # We can check if it's valid as-is, but the prompt asks for the value of 't'.
        # Since 't' doesn't exist, we conclude 'none'.
        print(f"Result for pattern {pattern_id} is 'none'.\n")
        return "none"

    # 3. Solve for 't' using angle sum = 360
    print("Step 2: Angle Sum Check (Sum of angles must be 360 degrees)")
    angles = []
    t_index = -1
    known_angle_sum = 0
    equation_str_parts = []
    for i, angle in enumerate(angles_str):
        if angle == 't':
            angles.append('t')
            t_index = i
            equation_str_parts.append('t')
        else:
            val = float(angle)
            angles.append(val)
            known_angle_sum += val
            equation_str_parts.append(str(int(val)))
            
    equation_str = " + ".join(equation_str_parts)
    t_value = 360.0 - known_angle_sum
    print(f"  Equation: {equation_str} = 360")
    print(f"  Solved for t: t = 360 - {known_angle_sum} = {t_value}")

    if t_value <= 0:
        print(f"  Check FAILED: Angle t must be positive, but t = {t_value}.")
        print(f"Result for pattern {pattern_id} is 'none'.\n")
        return "none"

    final_angles = angles[:]
    final_angles[t_index] = t_value

    # 4. Check Kawasaki's Theorem
    print("Step 3: Kawasaki's Theorem Check (Alternating angle sums must be 180)")
    s_odd = sum(final_angles[i] for i in range(0, len(final_angles), 2))
    s_even = sum(final_angles[i] for i in range(1, len(final_angles), 2))

    odd_indices_angles = [str(int(final_angles[i])) if i != t_index else f"t({int(t_value)})" for i in range(0, len(final_angles), 2)]
    even_indices_angles = [str(int(final_angles[i])) if i != t_index else f"t({int(t_value)})" for i in range(1, len(final_angles), 2)]
    
    print(f"  Odd sum: {' + '.join(odd_indices_angles)} = {s_odd}")
    print(f"  Even sum: {' + '.join(even_indices_angles)} = {s_even}")

    # Use a small tolerance for floating point comparisons
    if math.isclose(s_odd, 180.0) and math.isclose(s_even, 180.0):
        print("  Check PASSED: Both sums are 180.")
        result = int(t_value) if t_value == int(t_value) else t_value
        print(f"Result for pattern {pattern_id} is '{result}'.\n")
        return result
    else:
        print("  Check FAILED: Alternating sums are not equal to 180.")
        print(f"Result for pattern {pattern_id} is 'none'.\n")
        return "none"

def main():
    patterns = [
        [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
        [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
        [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
        [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
    ]
    
    results = []
    for i, p in enumerate(patterns):
        result = solve_crease_pattern(p, i + 1)
        results.append(str(result))

    final_output = f"[{','.join(results)}]"
    print("Final comma-separated list of results:")
    print(final_output)

if __name__ == "__main__":
    main()