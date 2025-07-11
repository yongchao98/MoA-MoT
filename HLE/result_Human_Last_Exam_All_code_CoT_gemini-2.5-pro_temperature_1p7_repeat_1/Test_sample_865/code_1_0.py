import math

def solve_crease_pattern(pattern_data, pattern_num):
    """
    Analyzes a single-vertex crease pattern to find the unknown angle 't'
    that allows for flat-folding, and explains the reasoning.

    Args:
        pattern_data (list): A list of angles (or 't') and fold types (M/V).
        pattern_num (int): The identifier for the pattern (e.g., 1, 2, 3, 4).

    Returns:
        A string representing the value of 't' or "none".
    """
    # 1. Parse the pattern and identify key properties
    angles_str = [str(a) for a in pattern_data[0::2]]
    folds = pattern_data[1::2]
    num_m = folds.count('M')
    num_v = folds.count('V')
    t_index = -1
    if 't' in angles_str:
        t_index = angles_str.index('t')
    
    print(f"--- Analysis for Pattern {pattern_num}: {pattern_data} ---")

    # 2. Check Maekawa's Theorem: |#M - #V| must be 2
    if abs(num_m - num_v) != 2:
        print(f"Result: Not flat-foldable.")
        print(f"Reason: Maekawa's Theorem is not satisfied. The number of mountain folds ({num_m}) and valley folds ({num_v}) must differ by 2, but |{num_m} - {num_v}| = {abs(num_m - num_v)}.")
        return "none"

    # Handle case with no variable 't' to solve for. This applies to Pattern 2.
    if t_index == -1:
        print(f"Result: Not flat-foldable.")
        print(f"Reason: The pattern is fixed (no unknown 't'), but it fails the flat-folding criteria.")
        return "none"

    # Maekawa's Theorem is satisfied, so we proceed.
    print(f"Info: Maekawa's Theorem holds (|{num_m} - {num_v}| = 2).")

    # 3. Formulate and solve two equations for 't'
    
    # Equation A: From the geometric constraint that angles sum to 360 degrees.
    known_angles_sum = sum(float(a) for a in angles_str if a != 't')
    t_from_sum = 360.0 - known_angles_sum
    
    # Equation B: From Kawasaki's Theorem that alternating angle sums are equal.
    odd_sum_known = sum(float(angles_str[i]) for i in range(0, len(angles_str), 2) if i != t_index)
    even_sum_known = sum(float(angles_str[i]) for i in range(1, len(angles_str), 2) if i != t_index)

    if t_index % 2 == 0:  # 't' is an odd-indexed angle (alpha_1, alpha_3...)
        t_from_kawasaki = even_sum_known - odd_sum_known
    else:  # 't' is an even-indexed angle (alpha_2, alpha_4...)
        t_from_kawasaki = odd_sum_known - even_sum_known
        
    print(f"Info: Angle sum requires t = 360.0 - {known_angles_sum} = {t_from_sum}")
    print(f"Info: Kawasaki's theorem requires t = {t_from_kawasaki}")
        
    # 4. Check for consistency and report the result.
    if math.isclose(t_from_sum, t_from_kawasaki) and t_from_sum > 0:
        solution = t_from_sum if not t_from_sum.is_integer() else int(t_from_sum)
        print(f"Result: A valid angle exists, t = {solution}")
        
        # Build and print the final equation as requested
        lhs_parts = []
        rhs_parts = []
        for i, part_str in enumerate(angles_str):
            current_part = str(solution) if i == t_index else part_str
            if i % 2 == 0:
                lhs_parts.append(current_part)
            else:
                rhs_parts.append(current_part)
        print(f"Final Equation: {' + '.join(lhs_parts)} = {' + '.join(rhs_parts)}")
        
        return str(solution)
    else:
        print(f"Result: Not flat-foldable.")
        print(f"Reason: The value required by angle sum (t={t_from_sum}) and Kawasaki's theorem (t={t_from_kawasaki}) are not consistent or t is not positive.")
        return "none"

def main():
    """Main function to run the analysis on all patterns."""
    patterns = [
        [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
        [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
        [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
        [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
    ]

    results = []
    for i, p in enumerate(patterns):
        result = solve_crease_pattern(p, i + 1)
        results.append(result)
        print("") # Add a newline for readability

    # Final answer in the required format
    print("--------------------------------------------------")
    print("Final answer in comma-separated list format:")
    print(f"[{','.join(results)}]")

if __name__ == '__main__':
    main()