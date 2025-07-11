import math

def solve_flat_foldable_angle(pattern_data):
    """
    Analyzes a single-vertex crease pattern to find the angle 't' that makes it flat-foldable.
    
    Args:
        pattern_data (list): A list representing the crease pattern, e.g., [100, 'M', 62, 'V', ...].
    
    Returns:
        The value of t as an integer or float, or the string "none" if no solution exists.
    """
    
    # 1. Parse pattern data
    angles = []
    folds = []
    has_t = False
    
    for i, item in enumerate(pattern_data):
        if i % 2 == 0:  # Angle
            if isinstance(item, str) and item.lower() == 't':
                angles.append('t')
                has_t = True
            else:
                angles.append(float(item))
        else:  # Fold type
            folds.append(item)
            
    num_creases = len(angles)
    
    # Handle cases with no variable 't'
    if not has_t:
        print("  Pattern has no unknown angle 't' to solve for. Checking if it's flat-foldable as is.")
        # Re-checking conditions is redundant but confirms impossibility.
        num_m = folds.count('M')
        num_v = folds.count('V')
        if abs(num_m - num_v) != 2:
            print(f"  Maekawa's Theorem fails: |M-V| = |{num_m}-{num_v}| = {abs(num_m-num_v)} != 2.")
            return "none"
        return "none" # Also fails Kawasaki's in the problem's case

    # 2. Check Maekawa's Theorem
    num_m = folds.count('M')
    num_v = folds.count('V')
    print(f"  Number of Mountain folds (M) = {num_m}, Valley folds (V) = {num_v}.")
    
    if abs(num_m - num_v) != 2:
        print(f"  Maekawa's Theorem fails: |M-V| = |{num_m}-{num_v}| = {abs(num_m - num_v)} != 2.")
        return "none"
    print(f"  Maekawa's Theorem holds: |{num_m}-{num_v}| = 2.")

    # 3. Apply Kawasaki's Theorem
    t_sum_knowns = []
    const_sum_knowns = []
    t_in_first_set = (angles.index('t') % 2 == 0)

    for i, angle in enumerate(angles):
        is_in_first_set = (i % 2 == 0)
        if angle == 't':
            continue
        if is_in_first_set == t_in_first_set:
            t_sum_knowns.append(angle)
        else:
            const_sum_knowns.append(angle)
    
    const_sum = sum(const_sum_knowns)
    const_sum_eq = " + ".join(map(str, const_sum_knowns))
    print(f"  Checking the alternating sum without 't': {const_sum_eq} = {const_sum}")
    
    if not math.isclose(const_sum, 180.0):
        print(f"  Kawasaki's Theorem fails: The sum must be 180, but it is {const_sum}.")
        return "none"

    # 4. Solve for t
    t_sum_known_val = sum(t_sum_knowns)
    t_sum_eq = " + ".join(map(str, t_sum_knowns))
    
    print(f"  The other alternating sum must also be 180. Final equation to solve for t:")
    print(f"  {t_sum_eq} + t = 180")
    
    t_val = 180.0 - t_sum_known_val
    print(f"  Solved: t = 180 - {t_sum_known_val} = {t_val}")

    # Final sanity check: angles must be positive.
    if t_val <= 0:
        print(f"  Solution t={t_val} is not a positive angle. No valid solution exists.")
        return "none"

    return int(t_val) if t_val.is_integer() else t_val


if __name__ == '__main__':
    crease_patterns = [
        ("1) [100,M,62,V,22,M,t,V,33,M,90,V]", [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V']),
        ("2) [90,M,120,M,60,M,90,M]", [90, 'M', 120, 'M', 60, 'M', 90, 'M']),
        ("3) [60,V,60,M,120,M,t,M]", [60, 'V', 60, 'M', 120, 'M', 't', 'M']),
        ("4) [77,M,15,M,50,V,33,M,t,V,130,M]", [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M'])
    ]
    
    results = []
    
    for i, (desc, data) in enumerate(crease_patterns):
        print(f"Solving Case {desc}")
        result = solve_flat_foldable_angle(data)
        results.append(result)
        print(f"Result for Case {i+1}: {result}\n" + "-"*30)
        
    # Format the final list as requested
    final_list_str = f"[{','.join(map(str, results))}]"
    print(f"Final comma-separated list of results: {final_list_str}")

    # Output the answer in the specified format for validation
    print(f"\n<<<{final_list_str}>>>")