import math

def solve_and_explain(pattern_list, case_number):
    """
    Analyzes a single-vertex crease pattern, prints the step-by-step reasoning,
    and returns the value of 't' or 'none'.
    """
    print(f"--- Analysis for Case {case_number}: {pattern_list} ---")

    # 1. Parse data from the input list
    angles_with_t = pattern_list[0::2]
    folds = pattern_list[1::2]
    num_creases = len(angles_with_t)

    # 2. Check for presence of 't'
    has_t = 't' in angles_with_t
    if not has_t:
        print("This pattern does not contain an unknown angle 't' to solve for.")
        m_count = folds.count('M')
        v_count = folds.count('V')
        if abs(m_count - v_count) != 2:
            print(f"The pattern also fails Maekawa's Theorem: |#M - #V| = |{m_count} - {v_count}| = {abs(m_count - v_count)}, which is not 2.")
        print("Result: none\n")
        return 'none'

    # 3. Check Maekawa's Theorem: |#M - #V| = 2
    m_count = folds.count('M')
    v_count = folds.count('V')
    maekawa_diff = abs(m_count - v_count)
    if maekawa_diff != 2:
        print(f"Failed Maekawa's Theorem: |#M - #V| = |{m_count} - {v_count}| = {maekawa_diff}, which is not 2.")
        print("Result: none\n")
        return 'none'
    else:
        print(f"Maekawa's Theorem holds: |#M - #V| = |{m_count} - {v_count}| = {maekawa_diff}.")

    # 4. Apply Kawasaki's Theorem to find 't'
    t_index = angles_with_t.index('t')
    known_angles = [float(a) for a in angles_with_t if a != 't']

    # 4a. From the total sum of angles being 360 degrees
    sum_of_knowns = sum(known_angles)
    t_from_sum = 360.0 - sum_of_knowns
    
    known_angles_str = [str(int(a)) for a in known_angles]
    sum_of_knowns_str = str(int(sum_of_knowns))
    t_from_sum_str = str(int(t_from_sum)) if t_from_sum.is_integer() else str(t_from_sum)

    print(f"Condition 1 (Total Sum = 360): {' + '.join(known_angles_str)} + t = 360")
    print(f"  => {sum_of_knowns_str} + t = 360")
    print(f"  => t = {t_from_sum_str}")

    # 4b. From the alternating sum of angles being equal
    sum_odd_pos_known = 0
    sum_even_pos_known = 0
    odd_pos_terms = []
    even_pos_terms = []

    for i in range(num_creases):
        term_val = angles_with_t[i]
        term_str = str(term_val) if isinstance(term_val, int) else term_val
        
        if i % 2 == 0:  # Odd positions (1st, 3rd, ...)
            odd_pos_terms.append(term_str)
            if term_val != 't':
                sum_odd_pos_known += float(term_val)
        else:  # Even positions (2nd, 4th, ...)
            even_pos_terms.append(term_str)
            if term_val != 't':
                sum_even_pos_known += float(term_val)
    
    t_from_alt = 0.0
    t_is_in_odd_pos = (t_index % 2 == 0)

    left_side = " + ".join(odd_pos_terms)
    right_side = " + ".join(even_pos_terms)
    print(f"Condition 2 (Alternating Sums Equal): {left_side} = {right_side}")
    
    sum_odd_pos_known_str = str(int(sum_odd_pos_known))
    sum_even_pos_known_str = str(int(sum_even_pos_known))
    
    if t_is_in_odd_pos:
        t_from_alt = sum_even_pos_known - sum_odd_pos_known
        print(f"  => {sum_odd_pos_known_str} + t = {sum_even_pos_known_str}")
    else:
        t_from_alt = sum_odd_pos_known - sum_even_pos_known
        print(f"  => {sum_odd_pos_known_str} = {sum_even_pos_known_str} + t")

    t_from_alt_str = str(int(t_from_alt)) if t_from_alt.is_integer() else str(t_from_alt)
    print(f"  => t = {t_from_alt_str}")

    # 5. Compare results and determine final answer
    if math.isclose(t_from_sum, t_from_alt) and t_from_sum > 0:
        result_val = int(t_from_sum)
        print(f"The values for t are consistent ({t_from_sum_str}) and positive.")
        print(f"Result: {result_val}\n")
        return result_val
    else:
        if not math.isclose(t_from_sum, t_from_alt):
            print(f"The values for t ({t_from_sum_str} and {t_from_alt_str}) from the two conditions are not consistent.")
        if t_from_sum <= 0:
            print(f"The calculated angle t={t_from_sum_str} is not a positive value.")
        print("Result: none\n")
        return 'none'

# --- Main Script Execution ---
patterns_data = [
    [100, 'M', 62, 'V', 22, 'M', 't', 'V', 33, 'M', 90, 'V'],
    [90, 'M', 120, 'M', 60, 'M', 90, 'M'],
    [60, 'V', 60, 'M', 120, 'M', 't', 'M'],
    [77, 'M', 15, 'M', 50, 'V', 33, 'M', 't', 'V', 130, 'M']
]

final_results = []
for i, pattern in enumerate(patterns_data):
    result = solve_and_explain(pattern, i + 1)
    final_results.append(result)

print("--- Final Answer ---")
print("The four values in order in a comma separated list are:")
final_result_str = ",".join(map(str, final_results))
print(f"[{final_result_str}]")
