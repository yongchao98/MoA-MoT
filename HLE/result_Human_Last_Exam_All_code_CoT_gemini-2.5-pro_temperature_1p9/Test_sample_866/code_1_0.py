import math

def solve_crease_pattern(data):
    """
    Calculates the number of valid flat-foldable assignments for a given 
    single-vertex crease pattern.
    """
    # 1. Parse Input
    angles = [x for x in data if isinstance(x, int)]
    creases = [x for x in data if isinstance(x, str)]

    n_creases = len(creases)
    print(f"Analyzing crease pattern with {n_creases} creases and {len(angles)} angles.")
    print("Angles:", angles)
    print("Crease assignments:", creases)
    print("-" * 30)

    # 2. Verify Geometric Constraints
    
    # Rule 1: Even number of creases
    # Maekawa's theorem implies N_M + N_V must be even.
    print("Rule 1: The total number of creases must be even.")
    is_even_crease_count = (n_creases % 2 == 0)
    if not is_even_crease_count:
        print(f"Check FAILED: The number of creases is {n_creases}, which is odd.")
        # We can stop here, but for completeness, we check the other rules.
    else:
        print(f"Check PASSED: The number of creases is {n_creases}, which is even.")

    # Rule 2: Maekawa's Theorem
    print("\nRule 2: Maekawa's Theorem (|N_M - N_V| = 2).")
    n_mountain_init = creases.count('M')
    n_valley_init = creases.count('V')
    n_unknown = creases.count('?')
    
    print(f"Initial counts: Mountains={n_mountain_init}, Valleys={n_valley_init}, Unknown={n_unknown}")
    print("Let 'm' be the number of '?' assigned as Mountain folds.")
    print(f"Then number of '?' assigned as Valley folds will be {n_unknown} - m.")
    print("The condition to satisfy is |(N_M_initial + m) - (N_V_initial + v)| = 2")
    
    # Deriving the equation for m
    # |(n_mountain_init + m) - (n_valley_init + n_unknown - m)| = 2
    # |n_mountain_init - n_valley_init - n_unknown + 2m| = 2
    c = n_mountain_init - n_valley_init - n_unknown
    
    print(f"\nThis gives the equation: |({n_mountain_init} - {n_valley_init} - {n_unknown}) + 2m| = 2, which is |{c} + 2m| = 2")

    # Case a: N_M - N_V = 2
    print(f"Checking case N_M - N_V = 2:  {c} + 2m = 2")
    rhs1 = 2 - c
    print(f"2m = {rhs1}")
    m1 = rhs1 / 2
    print(f"m = {m1}")

    # Case b: N_M - N_V = -2
    print(f"Checking case N_M - N_V = -2: {c} + 2m = -2")
    rhs2 = -2 - c
    print(f"2m = {rhs2}")
    m2 = rhs2 / 2
    print(f"m = {m2}")
    
    valid_m_solutions = []
    if m1.is_integer() and 0 <= m1 <= n_unknown:
        valid_m_solutions.append(int(m1))
    if m2.is_integer() and 0 <= m2 <= n_unknown and int(m2) not in valid_m_solutions:
        valid_m_solutions.append(int(m2))
    
    if not valid_m_solutions:
        print("\nCheck FAILED: No valid integer solution for 'm' exists.")
    else:
        print(f"\nCheck PASSED: Found valid integer solution(s) for m: {valid_m_solutions}")
    
    # Rule 3: Kawasaki's Theorem
    print("\nRule 3: Kawasaki's Theorem (Sum of alternating angles = 180 degrees).")
    kawasaki_ok = False
    if not is_even_crease_count:
        print("Check FAILED: The number of angles is odd, so Kawasaki's theorem is not applicable in its standard form.")
    else:
        sum_odd = sum(angles[0::2])
        sum_even = sum(angles[1::2])
        print(f"Sum of odd-indexed angles ({angles[0::2]}) = {sum_odd}")
        print(f"Sum of even-indexed angles ({angles[1::2]}) = {sum_even}")
        if sum_odd == 180 and sum_even == 180:
            print("Check PASSED: Angle sums are both 180.")
            kawasaki_ok = True
        else:
            print("Check FAILED: Angle sums are not both 180.")
    
    # 3. Final Conclusion and Calculation
    print("-" * 30)
    print("Final Conclusion:")
    if is_even_crease_count and kawasaki_ok and valid_m_solutions:
        total_assignments = 0
        for m in valid_m_solutions:
            # Number of ways to choose m creases to be 'M' out of n_unknown
            assignments = math.comb(n_unknown, m)
            print(f"For m={m}, number of assignments is C({n_unknown}, {m}) = {assignments}")
            total_assignments += assignments
    else:
        print("The pattern is not flat-foldable as it violates one or more fundamental rules.")
        total_assignments = 0

    print(f"\nTotal number of flat-foldable assignments: {total_assignments}")
    return total_assignments

# The input crease pattern from the user
input_pattern = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

# Run the solver and get the final answer
final_answer = solve_crease_pattern(input_pattern)