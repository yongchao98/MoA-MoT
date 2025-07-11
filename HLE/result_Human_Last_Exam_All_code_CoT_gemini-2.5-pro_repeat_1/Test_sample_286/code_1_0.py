import math

def solve_M():
    """
    Solves for the smallest positive integer M based on the problem's conditions.

    The problem reduces to solving the inequality W_N <= 2 * (1 - 10**(-M)),
    where W_N is the worst-case minimum range of partial sums, which is 1.
    """

    W_N = 1
    
    print("The problem requires finding the smallest positive integer M such that for any a_i in [0,1],")
    print("we can construct the sequence x_i in [-L, L] where L = 1 - 10**(-M).")
    print("This condition simplifies to the inequality:")
    print("sup_a inf_eps (max_k S_k - min_k S_k) <= 2 * L")
    print("\nThe value of the worst-case minimal range, sup_a inf_eps (max_k S_k - min_k S_k), is 1.")
    print("So, the inequality becomes:")
    
    # Using python's f-string formatting to show the equation with values
    print(f"{W_N} <= 2 * (1 - 10**(-M))")
    
    print("\nSolving for M:")
    print("1/2 <= 1 - 10**(-M)")
    print("10**(-M) <= 1 - 1/2")
    print("10**(-M) <= 1/2")
    print("1 / (10**M) <= 1 / 2")
    print("10**M >= 2")

    print("\nTaking the logarithm base 10 on both sides:")
    print("M >= log10(2)")

    log10_2 = math.log10(2)
    print(f"\nThe value of log10(2) is approximately: {log10_2}")

    # Find the smallest integer M
    # M must be a positive integer, and M >= 0.301...
    # The smallest such integer is 1.
    final_M = math.ceil(log10_2)
    if final_M == 0: # M must be a positive integer
        final_M = 1
        
    print(f"\nThe smallest positive integer M satisfying M >= {log10_2} is {final_M}.")
    
    # The final answer in the required format
    print(f"\nFinal Answer: {final_M}")

solve_M()