import math

def solve_exam_problem():
    """
    Solves the two-part combinatorial problem regarding exam creation.
    """
    # --- Part 1: n = 14, find maximum m (number of exams) ---
    print("--- Solving Part 1: If n = 14, how many different exams can be created? ---")
    print("Let n be the number of questions, k=4 be the number of questions per exam, and m be the number of exams.")
    print("The condition is that any two exams share at most one question: |E_i intersect E_j| <= 1.\n")

    print("Step 1: Find an upper bound for m using a pair-counting argument (Johnson bound).")
    n_part1 = 14
    k = 4
    pairs_per_exam = math.comb(k, 2)
    print(f"Each exam has C(k, 2) = C({k}, 2) = {pairs_per_exam} pairs of questions.")
    print("Since two exams can share at most one question, they cannot share any pair of questions.")
    print("Thus, all pairs of questions across all m exams must be distinct.")
    print("The total number of available pairs from n questions is C(n, 2).")
    total_pairs_part1 = math.comb(n_part1, 2)
    print(f"For n = {n_part1}, the total number of pairs is C({n_part1}, 2) = ({n_part1} * {n_part1-1}) / 2 = {total_pairs_part1}.")
    print(f"The total number of pairs from m exams is m * C(4, 2) = m * {pairs_per_exam}.")
    print(f"So, m * {pairs_per_exam} <= {total_pairs_part1}.")
    m_bound1 = math.floor(total_pairs_part1 / pairs_per_exam)
    print(f"m <= {total_pairs_part1} / {pairs_per_exam} => m <= {total_pairs_part1 / pairs_per_exam:.2f}")
    print(f"So, the number of exams m must be less than or equal to {m_bound1}.\n")

    print("Step 2: Refine the upper bound using the Schönheim bound.")
    print("The Schönheim bound states m <= floor(n/k * floor((n-1)/(k-1))).")
    v_minus_1_div_k_minus_1 = math.floor((n_part1 - 1) / (k - 1))
    print(f"For our values, floor(({n_part1}-1)/({k}-1)) = floor({n_part1-1}/{k-1}) = floor({(n_part1-1)/(k-1):.2f}) = {v_minus_1_div_k_minus_1}.")
    m_bound2 = math.floor(n_part1 / k * v_minus_1_div_k_minus_1)
    print(f"m <= floor({n_part1}/{k} * {v_minus_1_div_k_minus_1}) = floor({n_part1/k:.2f} * {v_minus_1_div_k_minus_1}) = floor({n_part1/k * v_minus_1_div_k_minus_1}) = {m_bound2}.")
    print(f"This gives a tighter bound: m <= {m_bound2}.\n")

    print("Step 3: Check the feasibility of the values at the bound.")
    print("For m = 14: If n=14 and m=14, the average number of times each question is used is (m*k)/n = (14*4)/14 = 4.")
    print("This implies every question must appear in exactly 4 exams. This structure, a (14,4,1)-configuration, is known in design theory to be impossible to construct.")
    print("Therefore, m cannot be 14.\n")

    print("For m = 13: We can show this is possible using a known combinatorial object called a projective plane of order 3, or PG(2,3).")
    print("This structure consists of exactly 13 points (questions) and 13 lines (exams).")
    print("In this structure, each exam contains 4 questions, and any two exams intersect at exactly one question, which satisfies our condition.")
    print("To use this for n=14, we can simply use questions {1, 2, ..., 13} from the pool of 14 and build the 13 exams according to this structure.")
    print("Thus, m = 13 is possible for n=14.\n")

    answer1 = 13
    print(f"Conclusion for Part 1: The maximum number of different exams is {answer1}.\n")
    print("-------------------------------------------\n")

    # --- Part 2: m = 10, find minimum n ---
    print("--- Solving Part 2: What is the minimum value of n needed to prepare 10 exams? ---")
    print("Here m = 10 exams and k = 4 questions per exam. We need to find the minimum n.\n")

    print("Step 1: Use a question frequency argument to find a lower bound for n.")
    m_part2 = 10
    total_slots = m_part2 * k
    print(f"The total number of question slots across all {m_part2} exams is m * k = {m_part2} * {k} = {total_slots}.")
    print("Let r_i be the number of exams question i appears in. The sum of all r_i over all n questions must be 40.")
    print("Consider any single question q. It appears in r_q exams. The other 3 questions in each of these r_q exams must come from the remaining n-1 questions, and their sets must be disjoint.")
    print(f"This leads to the inequality: r_q * (k-1) <= n-1, which is r_q * 3 <= n-1.\n")

    print("Step 2: Test if n=12 is possible.")
    n_test_12 = 12
    r_max_12 = math.floor((n_test_12 - 1) / (k - 1))
    print(f"If n = {n_test_12}, the maximum frequency for any question is r_max <= floor(({n_test_12}-1)/3) = {r_max_12}.")
    max_total_slots_12 = n_test_12 * r_max_12
    print(f"If every question appears at most {r_max_12} times, the total number of question slots is at most n * r_max = {n_test_12} * {r_max_12} = {max_total_slots_12}.")
    print(f"This contradicts the requirement that the total number of slots must be {total_slots}.")
    print(f"Since {max_total_slots_12} < {total_slots}, n = {n_test_12} is not possible.\n")

    print("Step 3: Test if n=13 is possible.")
    n_test_13 = 13
    print(f"If n = {n_test_13}, the max frequency is r_max <= floor(({n_test_13}-1)/3) = 4. The total slots could be up to 13 * 4 = 52, which is consistent with 40.")
    print("To show n=13 is sufficient, we can again use the projective plane of order 3 (PG(2,3)).")
    print("This design provides 13 exams on 13 questions that satisfy the conditions.")
    print(f"We only need to prepare 10 exams, so we can simply choose any 10 out of the 13 available exams from this design.")
    print(f"This provides a valid construction for m = 10 with n = {n_test_13}.")
    
    answer2 = 13
    print(f"\nConclusion for Part 2: Since n=12 is not possible, the minimum value of n is {answer2}.\n")

    print("================= Final Answer Summary =================")
    print(f"1. For n = 14, the maximum number of different exams is: {answer1}")
    print(f"2. To prepare 10 exams, the minimum n needed is: {answer2}")
    print("========================================================")

solve_exam_problem()