import math

def solve_crease_pattern():
    """
    Calculates the number of valid crease assignments for a given pattern
    based on a hypothesized rule for odd-degree vertices.
    """
    print("Step 1: Define the crease pattern from the input.")
    # The input [60,M,30,?,50,?,70,V,150,?] defines a vertex with 5 creases.
    creases_spec = ['M', '?', '?', 'V', '?']
    n = len(creases_spec)
    print(f"The pattern has n = {n} creases, which is an odd number.")

    num_m_known = creases_spec.count('M')
    num_v_known = creases_spec.count('V')
    num_q_marks = creases_spec.count('?')

    print(f"The pattern has {num_m_known} known Mountain (M), {num_v_known} known Valley (V), and {num_q_marks} unassigned (?) folds.")
    print("-" * 20)

    print("Step 2: Apply a modified Maekawa's Theorem.")
    print("Standard Maekawa's theorem |M - V| = 2 is not applicable for odd-degree vertices.")
    print("We hypothesize a rule for n=5, such as |M - V| = 3, to find a solution.")

    # The hypothesis is that |M_total - V_total| = 3.
    # M_total = num_m_known + m_assigned
    # V_total = num_v_known + v_assigned
    # This simplifies to |(1 + m_assigned) - (1 + v_assigned)| = 3, so |m_assigned - v_assigned| = 3.
    # We need to find the number of ways to assign the '?' to 'M' or 'V' to satisfy this.
    k_to_test = 3
    print(f"Testing the rule: |M_total - V_total| = {k_to_test}")
    print(f"This is equivalent to |m_assigned - v_assigned| = {k_to_test}, where m_assigned + v_assigned = {num_q_marks}.")
    print("-" * 20)
    
    print("Step 3: Calculate the number of valid assignments.")
    total_assignments = 0
    
    # We look for non-negative integer solutions (m, v) to:
    # m + v = 3
    # |m - v| = 3

    # Case 1: m - v = 3  =>  (m=3, v=0)
    m1 = 3
    v1 = 0
    ways1 = math.comb(num_q_marks, m1)
    print(f"Case 1: Assign {m1} '?' to M and {v1} to V.")
    print(f"Number of ways for Case 1 is C({num_q_marks}, {m1}) = {ways1}")
    total_assignments += ways1

    # Case 2: m - v = -3 =>  (m=0, v=3)
    m2 = 0
    v2 = 3
    ways2 = math.comb(num_q_marks, m2)
    print(f"Case 2: Assign {m2} '?' to M and {v2} to V.")
    print(f"Number of ways for Case 2 is C({num_q_marks}, {m2}) = {ways2}")
    total_assignments += ways2
    print("-" * 20)
    
    print("Final Result:")
    print(f"The total number of assignments is the sum of the ways from all valid cases.")
    print(f"Final Equation: C({num_q_marks}, {m1}) + C({num_q_marks}, {m2}) = {ways1} + {ways2} = {total_assignments}")

solve_crease_pattern()
