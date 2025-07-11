def print_complexity_analysis():
    """
    Prints a step-by-step analysis of the given algorithm's time complexity.
    """
    analysis_text = """
Step-by-step Analysis:

1.  **Algorithm Structure:**
    The provided algorithm is a recursive function `sort(a, n)`.
    -   Base Case: The recursion terminates when `n <= 1`. An array of size 1 or 0 is already considered sorted. This takes constant time, which we can denote as T(1) = c.
    -   Recursive Step: For `n > 1`, the function performs three main actions:
        1.  It calls itself with a smaller problem size: `sort(a, n - 1)`.
        2.  It performs a constant number of operations: one comparison (`if a[n - 1] < a[n - 2]`) and one potential swap (`swap(...)`).
        3.  It calls itself again with the same smaller problem size: `sort(a, n - 1)`.

2.  **Recurrence Relation Formulation:**
    Let T(n) represent the time complexity for `sort(a, n)`. Based on the structure above, we can write a recurrence relation. The total time is the sum of the times for its parts:
    -   T(n) = (Time for first recursive call) + (Time for comparison/swap) + (Time for second recursive call)
    -   T(n) = T(n - 1) + O(1) + T(n - 1)
    -   This simplifies to: T(n) = 2 * T(n - 1) + c (where c is a constant representing the O(1) work).

3.  **Solving the Recurrence Relation:**
    We can solve this by repeatedly expanding the T(n - 1) term:
    -   T(n) = 2 * T(n - 1) + c
    -   T(n) = 2 * (2 * T(n - 2) + c) + c = 4 * T(n - 2) + 2c + c
    -   T(n) = 4 * (2 * T(n - 3) + c) + 3c = 8 * T(n - 3) + 4c + 2c + c
    -   After i steps, a pattern emerges: T(n) = 2^i * T(n - i) + c * (2^(i-1) + ... + 2^1 + 2^0)

    The sum is a geometric series equal to 2^i - 1. So, T(n) = 2^i * T(n - i) + c * (2^i - 1).

    We continue until we hit the base case T(1), which happens when n - i = 1, or i = n - 1. Substituting i = n - 1:
    -   T(n) = 2^(n-1) * T(1) + c * (2^(n-1) - 1)

4.  **Determining the Big-O Complexity:**
    The final equation for T(n) is `2^(n-1) * T(1) + c * 2^(n-1) - c`.
    Since T(1) and c are constants, the term that grows fastest as n increases is 2^(n-1).
    In Big-O notation, we ignore constant coefficients and lower-order terms.
    Therefore, the time complexity is O(2^(n-1)), which is equivalent to O(2^n).

5.  **Conclusion:**
    The time complexity of this algorithm is O(2^n). This exponential growth makes the algorithm highly inefficient for anything but very small inputs.
"""
    print(analysis_text)

if __name__ == '__main__':
    print_complexity_analysis()