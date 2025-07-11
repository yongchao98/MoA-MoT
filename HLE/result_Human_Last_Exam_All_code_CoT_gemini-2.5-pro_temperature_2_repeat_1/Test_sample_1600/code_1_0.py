import math

def double_factorial(n):
    """Computes the double factorial n!!."""
    if n < 0:
        return 0
    if n == 0:
        return 1
    result = 1
    for i in range(n, 0, -2):
        result *= i
    return result

def solve():
    """
    Calculates a(3), the number of non-vanishing Feynman diagrams of order 6.

    a(n): number of irreducible pairings of 2n points on a line.
    A(n): total number of pairings of 2n points, given by (2n-1)!!.

    Recurrence relation: a(n) = A(n) - sum_{k=1 to n-1} a(k) * A(n-k)
    where A(0) is taken to be 1.
    """
    # Dictionary to store calculated values of a(n)
    a = {}
    # Dictionary to store calculated values of A(n)
    A = {0: 1}

    # Calculate for n=1
    n = 1
    A[n] = double_factorial(2 * n - 1)
    # a(1) = A(1)
    a[n] = A[n]

    # Calculate for n=2
    n = 2
    A[n] = double_factorial(2 * n - 1)
    # a(2) = A(2) - a(1)*A(1)
    a[n] = A[n] - (a[1] * A[1])

    # Calculate for n=3
    n = 3
    A[n] = double_factorial(2 * n - 1)
    # a(3) = A(3) - (a(1)*A(2) + a(2)*A(1))
    sum_term = a[1] * A[2] + a[2] * A[1]
    a[n] = A[n] - sum_term

    # Print the calculation for a(3)
    # a(3) = 15 - (1 * 3 + 2 * 1) = 10
    print(f"The number of diagrams a(3) is calculated as follows:")
    print(f"a(3) = A(3) - (a(1) * A(2) + a(2) * A(1))")
    print(f"a(3) = {A[3]} - ({a[1]} * {A[2]} + {a[2]} * {A[1]})")
    final_result = A[3] - (a[1] * A[2] + a[2] * A[1])
    print(f"a(3) = {final_result}")

solve()
<<<10>>>