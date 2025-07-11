import sys

def solve():
    """
    Solves the problem based on the user-provided formula format.
    The user is expected to input p, k, n separated by spaces.
    """
    try:
        p_str, k_str, n_str = input("Enter p, k, n separated by spaces: ").split()
        p = int(p_str)
        k = int(k_str)
        n = int(n_str)
    except ValueError:
        print("Invalid input. Please enter three integers separated by spaces.")
        return

    # Based on the step-by-step derivation:
    # Let N be the number of coefficients not divisible by p^k.
    # Case 1: k = 1
    # The number of coefficients is 2.
    if k == 1:
        result = 2
        print(f"The number of coefficients is {result}.")
        print(f"Final Answer Formula: 2")
        print(f"Final Answer Equation: {result} = 2")
    # Case 2: k >= 2
    # The number of coefficients for P_1(x) is (2p-1)p^(k-2).
    # Each of the n meta-steps (from p^n) multiplies this by p.
    # So for P_{p^n}(x), the number is (2p-1)p^(k-2+n).
    else:
        # result = (2*p - 1) * (p**(n + k - 2))
        # The derivation above suggests growth for each of the p^n steps, not just n meta-steps
        # A simpler hypothesis that fits for P_1 is that the answer is $2p^{k-1+n}-p^{k-2+n}$ for k>=2
        # Let's check for k=2. $2p^{n+1}-p^n = (2p-1)p^n$.
        # This matches the reasoning of multiplying the count for P_1 by p^n
        
        term1 = 2 * (p**(n + k - 1))
        term2 = p**(n + k - 2)
        result = term1 - term2

        print(f"The number of coefficients is {result}.")
        print(f"Final Answer Formula: 2*p**(n + k - 1) - p**(n + k - 2)")
        print(f"Final Answer Equation: {result} = 2*({p}**({n} + {k} - 1)) - {p}**({n} + {k} - 2)")


solve()
