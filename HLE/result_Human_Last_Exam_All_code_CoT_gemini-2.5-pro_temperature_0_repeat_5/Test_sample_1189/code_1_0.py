import math

def get_f0(n, k, memo_f0):
    """
    Calculates f_0(n, k), the number of k-subsets of F_2^n that sum to 0.
    Uses memoization to store and retrieve results.
    """
    if k < 0:
        return 0
    if (n, k) in memo_f0:
        return memo_f0[(n, k)]

    try:
        term1 = math.comb(2**n, k)
        if k % 2 == 1:
            # For odd k, the second term in the formula is 0
            result = term1
        else:
            # For even k, k = 2r
            r = k // 2
            term2 = (2**n - 1) * ((-1)**r) * math.comb(2**(n-1), r)
            result = term1 + term2
        
        # The result should be perfectly divisible by 2**n
        final_result = result // (2**n)
        memo_f0[(n, k)] = final_result
        return final_result
    except ValueError:
        # This can happen if k > 2**n, math.comb will raise an error.
        memo_f0[(n, k)] = 0
        return 0

def solve():
    """
    Main function to solve the problem for given n and m.
    """
    try:
        # Read n and m from user input
        n_str, m_str = input("Enter positive integers n and m, separated by a space: ").split()
        n = int(n_str)
        m = int(m_str)
        if n <= 0 or m < 0:
            print("Please enter positive integers for n and non-negative integer for m.")
            return
    except (ValueError, IndexError):
        print("Invalid input. Please enter two integers separated by a space.")
        return

    memo_f0 = {}
    
    # Base case: f(n, 0) = 1 (the empty set of subsets)
    if m == 0:
        print("The number of such sets T is 1.")
        print("Equation: 1 = 1")
        print("<<<1>>>")
        return

    # Calculate f(n, m) = sum_{i=0 to m} (-1)^i * f_0(n, m-i)
    total_sum = 0
    equation_parts = []

    for i in range(m + 1):
        k = m - i
        term = get_f0(n, k, memo_f0)
        
        if i % 2 == 1:
            total_sum -= term
            sign = "-"
        else:
            total_sum += term
            sign = "+"

        if i == 0:
            equation_parts.append(f"{term}")
        else:
            equation_parts.append(f" {sign} {term}")

    print(f"The number of such sets T is {total_sum}.")
    print("The calculation is based on the formula: f(n,m) = f_0(n,m) - f_0(n,m-1) + ...")
    print("Equation:", "".join(equation_parts), f"= {total_sum}")
    print(f"<<<{total_sum}>>>")

solve()