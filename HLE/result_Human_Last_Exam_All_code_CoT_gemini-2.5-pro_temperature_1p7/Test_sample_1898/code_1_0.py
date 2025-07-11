import math

def solve():
    """
    Solves the two parts of the math competition problem using an analytical approach.
    """
    k = 4

    print("--- Part 1: Finding the maximum number of exams for n=14 ---")
    n1 = 14
    print(f"Given n = {n1} questions and k = {k} questions per exam.\n")

    # Constraint 1: Based on counting pairs of questions
    # m * k * (k-1) <= n * (n-1)
    # m <= n * (n-1) / (k * (k-1))
    bound1 = math.floor(n1 * (n1 - 1) / (k * (k - 1)))
    print("Step 1: Apply the first constraint from counting question pairs.")
    print(f"The number of exams 'm' must satisfy: m * C(k,2) <= C(n,2)")
    print(f"m * C(4,2) <= C(14,2) ==> m * 6 <= 91")
    print(f"m <= 91 / 6 = {91/6:.2f}. Since m must be an integer, m <= {math.floor(91/6)}.")
    # A more general form is m <= n(n-1)/k(k-1)
    print(f"A more direct calculation: m <= n(n-1)/(k(k-1)) = 14*13/(4*3) = 182/12 = {182/12:.2f}")
    print(f"So, from this constraint, the maximum number of exams m <= {bound1}.\n")

    # Constraint 2: Based on question frequency
    # r_i <= (n-1)/(k-1) and m*k = sum(r_i) <= n * max(r_i)
    r_max = math.floor((n1 - 1) / (k - 1))
    bound2 = math.floor(n1 * r_max / k)
    print("Step 2: Apply the second constraint from question frequency.")
    print("Let r_i be the number of exams containing question i.")
    print(f"r_i must satisfy: r_i * (k-1) <= (n-1) ==> r_i * 3 <= 13 ==> r_i <= {13/3:.2f}.")
    print(f"Since r_i must be an integer, any question can appear in at most r_max = {r_max} exams.\n")
    
    print("Step 3: Combine constraints to find a tighter bound.")
    print("The total number of question slots across all exams is m*k.")
    print(f"This must equal the sum of frequencies: m*k = sum(r_i).")
    print(f"We know each r_i <= {r_max}, so sum(r_i) <= n * r_max.")
    print(f"This gives the inequality: m*k <= n * r_max")
    print(f"m * {k} <= {n1} * {r_max} ==> {k}m <= {n1*r_max} ==> m <= {n1*r_max/k}")
    print(f"So, the tighter bound is m <= {bound2}.\n")

    print("Step 4: Analysis of the m = 14 case.")
    print("If m = 14, the condition m*k <= n*r_max becomes an equality: 14*4 = 14*4 = 56.")
    print("This means every question must appear in exactly r = 4 exams.")
    print("This requires the existence of a specific combinatorial structure (a packing design D(14,4,2) of size 14).")
    print("It is a known result in advanced combinatorial mathematics that such a structure cannot be created.")
    print("Therefore, m = 14 is not possible.\n")

    print("Step 5: Conclusion for Part 1.")
    # For n=13, a design with m=13 exists (the projective plane of order 3).
    # So for n=14, we can certainly form 13 exams by ignoring one question.
    print("The maximum number of exams must be less than 14. The next integer is 13.")
    print("It is possible to construct 13 exams with 14 questions (a known design for n=13 exists, so it also exists for n=14).")
    ans1 = 13
    print(f"The maximum number of different exams is {ans1}.\n")

    print("--- Part 2: Finding the minimum number of questions for m=10 ---")
    m2 = 10
    print(f"Given m = {m2} exams and k = {k} questions per exam.\n")

    print("We need to find the smallest integer n that satisfies both constraints.")
    print("Constraint 1: n(n-1) >= m*k(k-1) = 10*4*3 = 120")
    print("Constraint 2: n * floor((n-1)/3) >= m*k = 10*4 = 40\n")
    
    n_min = k # We need at least k questions
    while True:
        # Check constraint 1
        cond1 = n_min * (n_min - 1) >= m2 * k * (k - 1)
        # Check constraint 2
        r_max_n = math.floor((n_min - 1) / (k - 1))
        cond2 = n_min * r_max_n >= m2 * k

        if not cond1:
            print(f"Testing n = {n_min}: Fails Constraint 1. ({n_min}*({n_min-1}) = {n_min*(n_min-1)} < 120)")
            n_min += 1
            continue
        
        if not cond2:
            print(f"Testing n = {n_min}: Constraint 1 holds ({n_min*(n_min-1)} >= 120).")
            print(f"Testing n = {n_min}: Fails Constraint 2. (n*floor((n-1)/3) = {n_min}*{r_max_n} = {n_min*r_max_n} < 40)")
            n_min += 1
            continue

        print(f"Testing n = {n_min}: Constraint 1 holds ({n_min*(n_min-1)} >= 120).")
        print(f"Testing n = {n_min}: Constraint 2 holds (n*floor((n-1)/3) = {n_min}*{r_max_n} = {n_min*r_max_n} >= 40).")
        print(f"Both constraints are satisfied for n = {n_min}.\n")
        break
    
    print("Step 6: Conclusion for Part 2.")
    print("The analytical bounds show n must be at least 13.")
    print("The existence of a design with n=13 questions and m=13 exams (projective plane of order 3) guarantees that we can construct 10 exams.")
    ans2 = n_min
    print(f"The minimum value of n needed is {ans2}.\n")

    print(f"Final Answer (Part 1, Part 2):")
    print(f"{ans1}, {ans2}")
    
solve()