import collections

def calculate_limit():
    """
    This function calculates the exact value of the limit lim_{N->inf} p(N)/N.
    
    The plan is as follows:
    1.  The expression for p(N) involves a sum over all possible coefficients (a, ..., g) and all possible values of m and n.
    2.  The key insight is to separate the coefficient tuples C=(a,...,g) into two cases.
        Case A: The polynomial P_C(x) = a*x^6 + ... + g is a constant function. This happens if and only if the coefficients of all non-constant terms are zero, i.e., a=b=c=d=e=f=0.
        Case B: The polynomial P_C(x) is not a constant function.
    3.  In Case A, the equation becomes F_n = -g, which is independent of m. For each of the N values of m (from 0 to N-1), we get the same number of solutions for n. The number of solutions for n, let's call it k_g, is constant for large N. Thus, each such coefficient tuple contributes N * k_g to p(N). The total contribution from Case A is N * (sum of k_g over g in [-25, 25]). This part of p(N) grows linearly with N.
    4.  In Case B, for a fixed coefficient tuple C, the equation F_n = -P_C(F_m) is a non-trivial Diophantine equation. It is known that such equations have only a finite number of integer solutions (n, m). Summing these finite numbers over all tuples in Case B gives a large but finite constant. This part of p(N) does not grow with N.
    5.  Therefore, p(N) is of the form K*N + C_1, where K = sum(k_g) and C_1 is a constant.
    6.  The limit lim_{N->inf} p(N)/N is then equal to K.
    7.  The code will calculate K by summing the number of solutions k_g for F_n = -g for each g from -25 to 25.
    """
    
    limit_g = 25

    # Generate Fibonacci numbers up to limit_g.
    fibs = [0, 1]
    while fibs[-1] <= limit_g:
        fibs.append(fibs[-1] + fibs[-2])

    # Count occurrences of each number in the Fibonacci sequence.
    # This directly gives the number of solutions for n for a given value.
    # F_1=1, F_2=1, so value 1 has 2 solutions (n=1, n=2).
    # Other Fibonacci numbers in this range are unique.
    fib_counts = collections.Counter(fibs)

    total_k = 0
    equation_terms = []

    print("The limit is determined by the cases where the polynomial is constant (a=b=c=d=e=f=0).")
    print("The equation simplifies to F_n = -g.")
    print("We sum the number of solutions for n for each g in [-25, 25]:\n")

    # Iterate through all possible integer values of g
    for g in range(-limit_g, limit_g + 1):
        value_fn = -g
        
        # F_n must be non-negative.
        if value_fn < 0:
            continue

        if value_fn in fib_counts:
            num_solutions_n = fib_counts[value_fn]
            if num_solutions_n > 0:
                print(f"g = {g:3d}  =>  F_n = {value_fn:2d}. Number of solutions for n: {num_solutions_n}")
                total_k += num_solutions_n
                equation_terms.append(str(num_solutions_n))

    print("\nThe limit is the sum of these solution counts.")
    # The order of summation doesn't matter, but sorting makes it canonical.
    equation_str = " + ".join(sorted(equation_terms, reverse=True))
    print(f"Final calculation: {equation_str} = {total_k}")
    
    print(f"\nThe exact value of the limit is {total_k}.")


calculate_limit()