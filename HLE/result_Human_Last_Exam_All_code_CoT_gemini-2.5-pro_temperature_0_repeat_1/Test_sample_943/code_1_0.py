import math

def solve_for_k():
    """
    This function calculates the value of K based on the derived equation.
    The slice rank of the tensor T_n was shown to be asymptotically 2^n.
    The problem provides the slice rank as (3 / 2^K)^n * e^(o(n)).
    By equating the asymptotic bases, we get the equation: 3 / (2^K) = 2.
    This script solves this equation for K.
    """

    # The numbers in the final equation
    val_3 = 3
    val_2a = 2
    val_2b = 2

    print("The equation to solve for K is derived by equating the asymptotic growth rates:")
    print(f"{val_3} / ({val_2a}^K) = {val_2b}")

    print("\nSolving the equation step-by-step:")
    print(f"1. Multiply both sides by {val_2a}^K: {val_3} = {val_2b} * {val_2a}^K")
    print(f"2. Combine the powers of 2: {val_3} = {val_2a}^(K + 1)")
    print(f"3. Take log base 2 of both sides: log2({val_3}) = K + 1")
    print(f"4. Isolate K: K = log2({val_3}) - 1")

    # Calculate the numerical value of K
    K = math.log2(val_3) - 1

    print(f"\nThe numerical value of K is:")
    print(K)

solve_for_k()