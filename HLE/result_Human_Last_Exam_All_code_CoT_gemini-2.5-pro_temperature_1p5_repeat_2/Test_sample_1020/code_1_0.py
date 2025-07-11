import math

def solve_complexity():
    """
    Calculates the complexity of a matrix product computing the polynomial f_n,d.

    The polynomial is f_{n,d} = sum_s x_{1,s(1)} * ... * x_{d,s(d)},
    where the sum is over all injective functions s: {1,...,d} -> {1,...,n}.

    A known construction for this polynomial uses a matrix product A_1 * ... * A_d
    where the intermediate dimensions m_i (for i=1 to d-1) are given by the
    number of subsets of size i from a set of n elements.
    - m_i = C(n, i) = n! / (i! * (n-i)!)

    The complexity of this matrix product is defined as 2 + m_1 + ... + m_{d-1}.
    This construction provides an upper bound for the minimal complexity. For certain
    values of n and d, it is known to be optimal or match other complex constructions.

    We will calculate the complexity for the example case n=7, d=3.
    """
    n = 7
    d = 3

    # Check if the conditions are met for the example
    if not (2 <= d <= n and d % 2 != 0):
        print(f"The example values n={n}, d={d} do not satisfy the problem constraints.")
        return

    # Calculate the complexity
    sum_of_m = 0
    for i in range(1, d):
        sum_of_m += math.comb(n, i)

    complexity = 2 + sum_of_m

    # Output the explanation and calculation
    print(f"For n = {n} and d = {d}, we calculate the complexity of computing the polynomial f_({n},{d}).")
    print("A general construction gives the complexity as C = 2 + sum_{i=1}^{d-1} C(n, i).")
    print("\nCalculation:")
    
    # Build the equation string with numbers
    equation_str = f"C = 2"
    for i in range(1, d):
        equation_str += f" + C({n}, {i})"
    
    print(equation_str)

    # Build the values string
    values_str = f"C = 2"
    for i in range(1, d):
        values_str += f" + {math.comb(n, i)}"
    
    print(values_str)
    
    # Print the final result
    print(f"C = {complexity}")

if __name__ == "__main__":
    solve_complexity()