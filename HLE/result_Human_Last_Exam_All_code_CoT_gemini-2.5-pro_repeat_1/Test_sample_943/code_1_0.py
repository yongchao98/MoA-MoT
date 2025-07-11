import math

def solve_for_k():
    """
    This function solves for the value of K based on the analysis.
    The slice rank of the tensor is 2^n - 1.
    The problem provides the formula (3/2^K)^n * e^(o(n)).
    Equating the asymptotic bases, we get the equation 2 = 3 / (2^K).
    This script solves this equation for K.
    """
    # The equation is 2 = 3 / (2^K)
    # which simplifies to 2^(K+1) = 3
    # K+1 = log2(3)
    # K = log2(3) - 1
    k_value = math.log2(3) - 1

    print("The analysis of the tensor's slice rank leads to the equation for the asymptotic base:")
    # As requested, printing the numbers in the final equation.
    print("2 = 3 / (2^K)")
    print("\nSolving for K:")
    print("2 * 2^K = 3")
    print("2^(K+1) = 3")
    print("K+1 = log2(3)")
    print("K = log2(3) - 1")
    print(f"\nThe numerical value of K is: {k_value}")

solve_for_k()
