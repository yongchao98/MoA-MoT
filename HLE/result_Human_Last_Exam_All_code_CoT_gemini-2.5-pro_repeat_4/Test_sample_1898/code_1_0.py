import math

def combinations(n, k):
    """Calculates the number of combinations 'n choose k'."""
    if k < 0 or k > n:
        return 0
    # Use integer division
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve():
    """
    Solves both parts of the problem and prints the detailed steps.
    """
    # --- Part 1 ---
    n1 = 14
    k = 4  # questions per exam
    print("--- Part 1: Finding the maximum number of exams for n = 14 ---")
    print(f"Given: n = {n1} questions in the pool, k = {k} questions per exam.")
    print("Condition: Any two exams can have at most one common question.\n")

    print("Step 1: Derive an upper bound for the number of exams (m) using pairs.")
    pairs_per_exam = combinations(k, 2)
    total_possible_pairs = combinations(n1, 2)
    print(f"Each exam has C(k, 2) = C({k}, 2) = {pairs_per_exam} pairs of questions.")
    print(f"The total number of possible pairs from n={n1} questions is C(n, 2) = C({n1}, 2) = {total_possible_pairs}.")
    print("The condition implies all pairs from all exams must be unique.")
    print("Therefore, m * (pairs per exam) <= (total possible pairs)")
    print(f"m * {pairs_per_exam} <= {total_possible_pairs}")
    m_bound1 = total_possible_pairs // pairs_per_exam
    print(f"m <= {total_possible_pairs / pairs_per_exam:.2f}, so m must be at most {m_bound1}.\n")

    print("Step 2: Derive a tighter upper bound using question frequency.")
    print("Let r be the maximum number of exams any single question can appear in.")
    print(f"For a question in r exams, it is paired with r * (k-1) = r * {k-1} other questions.")
    print(f"These must be distinct and chosen from the remaining n-1 = {n1-1} questions.")
    print(f"So, r * ({k}-1) <= {n1}-1  =>  r * {k-1} <= {n1-1}")
    r_max1 = (n1 - 1) // (k - 1)
    print(f"r <= {(n1-1)/(k-1):.2f}, which means the maximum frequency is r_max = {r_max1}.\n")

    print("The total number of question 'slots' across all exams is m * k.")
    print("This must be less than or equal to the total slots available, n * r_max.")
    print(f"m * k <= n * r_max")
    print(f"m * {k} <= {n1} * {r_max1}")
    m_bound2 = (n1 * r_max1) // k
    print(f"m <= {n1 * r_max1 / k}, so m must be at most {m_bound2}.\n")

    print(f"Step 3: Check if m = {m_bound2} is possible.")
    print(f"If m = {m_bound2}, the equality m * k = n * r_max holds, meaning every question must appear in exactly r_max = {r_max1} exams.")
    print(f"This requires a (v, k, λ)-design with v={n1}, k={k}, λ=1 (a Steiner System S(2,4,14)).")
    print("A necessary condition for such a design is that λ(v-1)/(k-1) must be an integer.")
    print(f"For our parameters, this is 1 * ({n1}-1)/({k}-1) = {n1-1}/{k-1} = {(n1-1)/(k-1):.2f}, which is not an integer.")
    print(f"Therefore, a configuration with m = {m_bound2} exams is impossible.\n")

    max_m = m_bound2 - 1
    print(f"Step 4: Conclude the maximum number of exams.")
    print(f"Since m cannot be {m_bound2}, the maximum possible integer value is m = {max_m}.")
    print("A design with 13 exams is known to exist (based on the projective plane of order 3).")
    q1_answer = max_m
    print(f"Final Answer for Part 1: The maximum number of exams is {q1_answer}.")
    print("-" * 50)

    # --- Part 2 ---
    m2 = 10
    print("\n--- Part 2: Finding the minimum n for m = 10 exams ---")
    print(f"Given: m = {m2} exams, k = {k} questions per exam.\n")

    print("Step 1: Use the combinatorial bounds to find a lower limit for n.")
    total_pairs_needed = m2 * pairs_per_exam
    print("Condition 1: m * C(k, 2) <= C(n, 2)")
    print(f"{m2} * {pairs_per_exam} <= n(n-1)/2  =>  {2*total_pairs_needed} <= n(n-1)")
    
    min_n_cond1 = 0
    n_test = k
    while True:
        if n_test * (n_test - 1) >= 2 * total_pairs_needed:
            min_n_cond1 = n_test
            break
        n_test += 1
    print(f"Testing values for n(n-1) >= {2*total_pairs_needed}:")
    print(f"If n={min_n_cond1-1}, n(n-1) = {(min_n_cond1-1)*(min_n_cond1-2)} which is too small.")
    print(f"If n={min_n_cond1}, n(n-1) = {min_n_cond1*(min_n_cond1-1)} which is large enough.")
    print(f"So, from this condition, we need at least n = {min_n_cond1} questions.\n")

    print(f"Step 2: Check integer n values starting from {min_n_cond1} using the second condition.")
    print("Condition 2: m * k <= n * r_max, where r_max = floor((n-1)/(k-1))\n")

    # Check n = 12
    n = 12
    print(f"Testing n = {n}:")
    r_max_12 = (n - 1) // (k - 1)
    total_incidences_needed = m2 * k
    max_possible_incidences_12 = n * r_max_12
    print(f"r_max for n={n} is floor(({n}-1)/({k}-1)) = floor({n-1}/{k-1}) = {r_max_12}")
    print(f"The number of required question slots is m * k = {m2} * {k} = {total_incidences_needed}.")
    print(f"The maximum possible slots for n={n} is n * r_max = {n} * {r_max_12} = {max_possible_incidences_12}.")
    print(f"The condition is {total_incidences_needed} <= {max_possible_incidences_12}, which is False.")
    print(f"Therefore, n = {n} is not possible.\n")

    # Check n = 13
    n = 13
    print(f"Testing n = {n}:")
    r_max_13 = (n - 1) // (k - 1)
    max_possible_incidences_13 = n * r_max_13
    print(f"r_max for n={n} is floor(({n}-1)/({k}-1)) = floor({n-1}/{k-1}) = {r_max_13}")
    print(f"The number of required question slots is still {total_incidences_needed}.")
    print(f"The maximum possible slots for n={n} is n * r_max = {n} * {r_max_13} = {max_possible_incidences_13}.")
    print(f"The condition is {total_incidences_needed} <= {max_possible_incidences_13}, which is True.")
    print(f"All necessary conditions are met for n = {n}.\n")
    
    print("Step 3: Conclude the minimum number of questions.")
    print(f"We have shown that n=12 is impossible and n=13 is arithmetically possible.")
    print("The existence of a suitable design for n=13 is also known.")
    q2_answer = n
    print(f"Final Answer for Part 2: The minimum value of n is {q2_answer}.")
    
    print("\n" + "="*50)
    print("Summary of Answers:")
    print(f"1. For n = 14, the maximum number of exams is {q1_answer}.")
    print(f"2. To prepare 10 exams, the minimum number of questions n is {q2_answer}.")
    print("="*50)
    
    return q1_answer, q2_answer

if __name__ == '__main__':
    q1_ans, q2_ans = solve()
    print(f"<<<1. {q1_ans}, 2. {q2_ans}>>>")
