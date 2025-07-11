def explain_complexity_analysis():
    """
    This function prints a step-by-step analysis of the provided
    algorithm's time complexity.
    """

    explanation = """
Let's analyze the time complexity of the given `sort(a, n)` algorithm.
Let T(n) be the time taken to run `sort` on an input of size `n`.

The pseudocode for the algorithm is:
def sort(a, n):
  if n > 1:
    sort(a, n - 1)
    if a[n - 1] < a[n - 2]:
      swap(a[n - 1], a[n - 2])
      sort(a, n - 1)

1.  **Base Case:** If n <= 1, the function does nothing and returns. This takes constant time, O(1).
    So, T(1) = c' (a constant).

2.  **Recursive Step (n > 1):** We analyze the operations:
    - `sort(a, n - 1)`: A recursive call on a subproblem of size n-1. Time taken: T(n-1).
    - `if a[n-1] < a[n-2]`: A comparison. Time taken: O(1).
    - `swap(...)`: A swap. Time taken: O(1).
    - `sort(a, n - 1)`: A second recursive call on a subproblem of size n-1. Time taken: T(n-1).

3.  **Worst-Case Scenario:** The worst case occurs when the condition `a[n - 1] < a[n - 2]` is always true, forcing the `swap` and the second recursive call to be executed at every step of the main recursion. This happens, for example, when the input array is sorted in reverse order.

4.  **Recurrence Relation:** In the worst case, the total time T(n) is the sum of the time for the two recursive calls plus the constant time for the comparison and swap.
    This gives us the recurrence relation:
    T(n) = T(n - 1) + T(n - 1) + c  (where c represents constant time work)
    T(n) = 2 * T(n - 1) + c

5.  **Solving the Recurrence Relation:** Let's solve this by expansion (unrolling the recursion).

    T(n) = 2 * T(n - 1) + c
    
    Substituting T(n-1) = 2 * T(n-2) + c:
    T(n) = 2 * (2 * T(n - 2) + c) + c
         = 4 * T(n - 2) + 2*c + c
         = 4 * T(n - 2) + 3*c

    Substituting T(n-2) = 2 * T(n-3) + c:
    T(n) = 4 * (2 * T(n - 3) + c) + 3*c
         = 8 * T(n - 3) + 4*c + 3*c
         = 8 * T(n - 3) + 7*c
    
    We can see a pattern here:
    T(n) = 2^k * T(n - k) + (2^k - 1) * c

    To find the final expression, we unroll until we hit the base case T(1). Let n - k = 1, so k = n - 1.
    T(n) = 2^(n-1) * T(1) + (2^(n-1) - 1) * c

    Since T(1) is a constant, the expression is dominated by the 2^(n-1) term.
    
    This shows that the runtime grows exponentially with `n`.
    Therefore, the time complexity is O(2^n).
"""
    print(explanation)

if __name__ == "__main__":
    explain_complexity_analysis()