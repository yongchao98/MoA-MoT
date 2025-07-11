import math

def solve():
    """
    Solves the two parts of the math competition exam problem.
    """

    # --- Part 1: n = 14, find max m ---
    print("--- Part 1: n = 14, find max m ---")
    n1 = 14
    k = 4
    
    print("\nStep 1: Analytical approach to find an upper bound for the number of exams (m).")
    print("Let m be the number of exams and n be the number of questions.")
    print("Each exam has k questions. Any two exams have at most 1 common question.")
    
    print("\nFirst, we consider pairs of questions.")
    print("Each exam consists of k=4 questions, so it contains C(k, 2) pairs of questions.")
    c_k_2 = math.comb(k, 2)
    print(f"C(4, 2) = (4 * 3) / 2 = {c_k_2}")
    
    print("\nAcross m exams, there are m * C(4, 2) such pairs.")
    print("Since any two exams can share at most one question, they cannot share a pair of questions.")
    print("Therefore, all these pairs must be distinct.")
    print("The total number of distinct pairs available from n=14 questions is C(n, 2).")
    c_n1_2 = math.comb(n1, 2)
    print(f"C(14, 2) = (14 * 13) / 2 = {c_n1_2}")
    
    print("\nThis leads to the inequality: m * C(4, 2) <= C(14, 2)")
    m_bound1 = c_n1_2 / c_k_2
    print(f"m * {c_k_2} <= {c_n1_2}  =>  m <= {c_n1_2 / c_k_2:.2f}")
    m_upper_bound1 = math.floor(m_bound1)
    print(f"So, m must be at most {m_upper_bound1}.")

    print("\nNext, we can find a tighter bound by considering individual questions.")
    print("Let r_j be the number of exams containing question j.")
    print("For any given question j, the r_j exams that contain it must not share any other question.")
    print("This means the sets of (k-1) other questions in each of these r_j exams are disjoint.")
    print("These disjoint sets are drawn from the remaining (n-1) questions.")
    print("Thus, r_j * (k-1) <= n-1.")
    r_max = math.floor((n1 - 1) / (k - 1))
    print(f"r_j * ({k}-1) <= {n1}-1  =>  r_j * {k-1} <= {n1-1}  => r_j <= {(n1-1)/(k-1):.2f}")
    print(f"This means any question can appear in at most r_max = {r_max} exams.")
    
    print("\nThe total number of question slots in all exams is m*k.")
    print("This is also the sum of r_j over all n questions.")
    print("m * k = sum(r_j) <= n * r_max")
    m_bound2 = (n1 * r_max) / k
    print(f"m * {k} <= {n1} * {r_max}  =>  m <= {n1 * r_max / k}")
    m_upper_bound2 = math.floor(m_bound2)
    print(f"This gives a tighter upper bound: m <= {m_upper_bound2}.")

    print("\nStep 2: Checking if m = 14 is achievable.")
    print("If m = 14, then every question must appear in exactly r_max = 4 exams.")
    print("This would form a (2,4,14)-design, also known as a Steiner System S(2,4,14).")
    print("A necessary condition for the existence of an S(2,k,v) is that v*(v-1) must be divisible by k*(k-1).")
    v, K = n1, k
    numerator = v * (v - 1)
    denominator = K * (K - 1)
    print(f"Check: Is {v}*({v}-1) divisible by {K}*({K}-1)?")
    print(f"({numerator} / {denominator}) = {numerator / denominator:.2f}, which is not an integer.")
    print("Since the condition is not met, an S(2,4,14) does not exist. Therefore, m cannot be 14.")
    
    print("\nStep 3: Concluding the maximum m.")
    print("From the analysis, m must be less than 14, so m <= 13.")
    print("To show that m = 13 is achievable, we note that a Steiner system S(2,4,13) does exist (the projective plane of order 3).")
    print("This system defines 13 exams on 13 questions satisfying the conditions.")
    print("With n=14 questions available, we can easily form these 13 exams using only 13 of the questions.")
    
    ans1 = 13
    print(f"\nConclusion for Part 1: The maximum number of different exams is {ans1}.")
    
    print("\n" + "="*50 + "\n")
    
    # --- Part 2: m = 10, find min n ---
    print("--- Part 2: m = 10, find min n ---")
    m2 = 10
    
    print("\nStep 1: Analytical approach to find a lower bound for the number of questions (n).")
    print("We use the same inequalities as in Part 1, but solve for n.")
    
    print("\nFrom the pairs inequality: m * C(k, 2) <= C(n, 2)")
    min_pairs = m2 * c_k_2
    print(f"{m2} * {c_k_2} <= n * (n - 1) / 2")
    print(f"{min_pairs} <= n * (n - 1) / 2")
    print(f"{2 * min_pairs} <= n * (n - 1)")
    
    n_lower_bound1 = 0
    while n_lower_bound1 * (n_lower_bound1 - 1) < 2 * min_pairs:
        n_lower_bound1 += 1
    print(f"By testing values (e.g., 11*10=110, 12*11=132), we find that n must be at least {n_lower_bound1}.")

    print("\nFrom the individual question inequality: m * k <= n * floor((n-1)/(k-1))")
    mk_val = m2 * k
    print(f"{m2} * {k} <= n * floor((n - 1) / ({k} - 1))")
    print(f"{mk_val} <= n * floor((n - 1) / 3)")
    
    print(f"We test values of n starting from our bound {n_lower_bound1}:")
    n_test = n_lower_bound1
    while True:
        rhs = n_test * math.floor((n_test - 1) / (k - 1))
        print(f"For n = {n_test}: n * floor((n-1)/3) = {n_test} * {math.floor((n_test-1)/3)} = {rhs}")
        if mk_val <= rhs:
            print(f"Since {mk_val} <= {rhs}, this condition is met. The minimum n is likely {n_test}.")
            n_lower_bound2 = n_test
            break
        else:
            print(f"Since {mk_val} > {rhs}, this condition is not met.")
            n_test += 1
    
    print(f"\nStep 2: Concluding the minimum n.")
    print(f"The necessary conditions require n to be at least {n_lower_bound2}.")
    print(f"To show n = {n_lower_bound2} is sufficient, we must show that 10 exams can be constructed.")
    print(f"As established before, the Steiner system S(2, 4, 13) exists, providing 13 exams on 13 questions.")
    print("Since we only need to prepare 10 exams, we can select 10 of these 13, and the condition will hold.")
    
    ans2 = 13
    print(f"\nConclusion for Part 2: The minimum value of n needed is {ans2}.")


if __name__ == '__main__':
    solve()