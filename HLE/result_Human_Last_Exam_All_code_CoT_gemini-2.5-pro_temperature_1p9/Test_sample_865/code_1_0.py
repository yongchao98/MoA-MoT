import re

def solve_flat_fold_angle(patterns):
    """
    Analyzes a list of single-vertex crease patterns to find the angle 't'
    that makes them flat-foldable.

    The function applies three main conditions for flat-foldability:
    1. Maekawa-Jardin Theorem: |#Mountain - #Valley folds| = 2
    2. Sum of all angles at the vertex must be 360 degrees.
    3. Kawasaki's Theorem: The sum of alternating angles must equal 180 degrees.

    It prints the step-by-step analysis for each pattern, including the
    equations used.
    """
    results = []

    for i, pattern_str in enumerate(patterns):
        print(f"--- Analyzing Case {i+1}: {pattern_str} ---")

        # 1. Parse the input string to get angles and folds
        parts = pattern_str.strip('[]').split(',')
        angles_raw = []
        folds = []
        t_index = -1
        is_t_present = False
        
        # The first item is an angle, so we iterate through pairs of (angle, fold)
        # after the first one.
        angles_raw.append(parts[0])
        for j in range(1, len(parts), 2):
            folds.append(parts[j])
            if (j + 1) < len(parts):
                angles_raw.append(parts[j+1])
        
        angles = []
        for idx, ang in enumerate(angles_raw):
            if ang == 't':
                angles.append('t')
                is_t_present = True
                t_index = idx
            else:
                angles.append(int(ang))

        # 2. Check Maekawa-Jardin Theorem
        m_count = folds.count('M')
        v_count = folds.count('V')
        if abs(m_count - v_count) != 2:
            print(f"Result: Fails Maekawa-Jardin Theorem.")
            print(f"Explanation: The difference between mountain ({m_count}) and valley ({v_count}) folds must be 2, but it is {abs(m_count - v_count)}.")
            results.append("none")
            print("-" * 20)
            continue
        
        # Handle cases with no 't'
        if not is_t_present:
            print(f"Result: No unknown angle 't' to solve for.")
            results.append("none")
            print("-" * 20)
            continue

        # 3. Use Sum=360 to find a candidate for 't'
        known_angles = [a for a in angles if a != 't']
        sum_known = sum(known_angles)
        
        t_candidate = 360 - sum_known
        
        equation_parts = [str(a) for a in known_angles]
        print(f"Step 1: Use the angle sum rule (Sum = 360).")
        print(f"Equation: {' + '.join(equation_parts)} + t = 360")
        print(f"Solving: {sum_known} + t = 360  =>  t = {t_candidate}")

        if t_candidate <= 0:
            print(f"Result: The calculated angle t={t_candidate} is not a valid positive angle.")
            results.append("none")
            print("-" * 20)
            continue

        # 4. Verify the candidate 't' with Kawasaki's Theorem
        final_angles = list(angles)
        final_angles[t_index] = t_candidate
        
        odd_angles = [final_angles[k] for k in range(0, len(final_angles), 2)]
        even_angles = [final_angles[k] for k in range(1, len(final_angles), 2)]
        
        odd_sum = sum(odd_angles)
        even_sum = sum(even_angles)
        
        odd_eq = ' + '.join(map(str, odd_angles))
        even_eq = ' + '.join(map(str, even_angles))

        print(f"Step 2: Verify with Kawasaki's Theorem using t = {t_candidate}.")
        print(f"Alternating Sum 1 (odd indices): {odd_eq} = {odd_sum}")
        print(f"Alternating Sum 2 (even indices): {even_eq} = {even_sum}")

        if odd_sum == 180 and even_sum == 180:
            print(f"Result: Success! t = {t_candidate} satisfies all conditions.")
            results.append(str(t_candidate))
        else:
            print(f"Result: Fails Kawasaki's Theorem. The alternating sums are not both 180.")
            print("Explanation: No value of 't' can simultaneously satisfy both the angle sum rule and Kawasaki's theorem.")
            results.append("none")
        print("-" * 20)
    
    # Print the final list
    print(f"\nFinal compiled list: [{','.join(results)}]")


if __name__ == '__main__':
    patterns_to_solve = [
        "[100,M,62,V,22,M,t,V,33,M,90,V]",
        "[90,M,120,M,60,M,90,M]",
        "[60,V,60,M,120,M,t,M]",
        "[77,M,15,M,50,V,33,M,t,V,130,M]"
    ]
    solve_flat_fold_angle(patterns_to_solve)