import math
import sys

def solve():
    """
    Calculates the smallest known complexity for computing the polynomial f_n,d.

    The problem asks for the smallest complexity of a matrix product computing f_{n,d},
    defined as the sum of x_{1,s(1)}*...*x_{d,s(d)} over all injective functions 
    s: {1,...,d} -> {1,...,n}. The complexity is 2 + m_1 + ... + m_{d-1}.

    A standard construction, conjectured to be optimal, uses intermediate dimensions
    m_i = C(n, i), the number of subsets of indices of size i. This is necessary to
    keep track of the indices used at each stage of the computation.

    This leads to the complexity formula: C = 2 + sum_{i=1}^{d-1} C(n, i).
    The provided code calculates this value.
    """
    
    # Check for command-line arguments, otherwise use default values for demonstration.
    if len(sys.argv) == 3:
        n_str = sys.argv[1]
        d_str = sys.argv[2]
        print(f"Using provided values: n={n_str}, d={d_str}\n")
    else:
        n_str = "7"
        d_str = "5"
        print("Usage: python this_script.py <n> <d>")
        print(f"Using default values for demonstration: n={n_str}, d={d_str}\n")

    try:
        n = int(n_str)
        d = int(d_str)

        if not (2 <= d <= n):
            print(f"Error: Constraints not met. Required: 2 <= d <= n, but got n={n}, d={d}.")
            return
        if d % 2 == 0:
            print(f"Error: Constraint not met. Required: d is odd, but got d={d}.")
            return

    except (ValueError, IndexError):
        print("Error: Invalid input. Please provide two integers for n and d.")
        return

    complexity_terms = []
    total_sum = 0
    for i in range(1, d):
        term = math.comb(n, i)
        complexity_terms.append(term)
        total_sum += term
        
    final_complexity = 2 + total_sum
    
    print("The smallest known complexity for computing f_n,d is given by the formula:")
    print("C(n,d) = 2 + sum_{i=1}^{d-1} C(n,i)\n")
    print(f"For n={n} and d={d}, the calculation is:")
    
    # We construct the equation string with each term listed explicitly.
    equation_parts = ["2"] + [str(t) for t in complexity_terms]
    equation = "Complexity = " + " + ".join(equation_parts)
    equation += f" = {final_complexity}"
    
    print(equation)
    
    print(f"\n<<< {final_complexity} >>>")

if __name__ == '__main__':
    solve()