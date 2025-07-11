def main():
    """
    Analyzes and explains the time complexity of the given recursive sorting algorithm.
    """
    print("### Time Complexity Analysis ###")
    
    print("\nThe given sorting algorithm is recursive. Let's analyze its structure:")
    print("def sort(a, n):")
    print("  if n > 1:")
    print("    sort(a, n - 1)          # Recursive call 1")
    print("    if a[n - 1] < a[n - 2]:  # Constant time operation")
    print("      swap(a[n-1], a[n - 2]) # Constant time operation")
    print("      sort(a, n - 1)        # Recursive call 2 (conditional)")
    
    print("\nLet T(n) be the time complexity for an input of size n.")
    
    print("\n--- Worst-Case Scenario ---")
    print("The worst case happens when the `if` condition is met at every step of the recursion.")
    print("This causes two recursive calls on a subproblem of size (n-1) to be executed.")
    print("This scenario leads to the following recurrence relation:")
    print("T(n) = T(n - 1) + T(n - 1) + O(1)")
    print("T(n) = 2 * T(n - 1) + c, where 'c' is a constant for the comparison and swap operations.")
    
    print("\n--- Solving the Recurrence Relation ---")
    print("We can solve this recurrence relation by repeatedly substituting the formula for T.")
    
    # Initial recurrence
    t1_term = 2
    c1_val = 1
    print("\nThe equation is:")
    print(f"  T(n) = {t1_term} * T(n - 1) + {c1_val} * c")
    
    # First expansion
    t2_term = t1_term * 2
    c2_val = t1_term * c1_val + c1_val
    print("\nAfter substituting T(n-1) = 2*T(n-2) + c:")
    print(f"  T(n) = {t1_term} * (2*T(n-2) + c) + {c1_val}*c")
    print(f"  T(n) = {t2_term} * T(n - 2) + {c2_val} * c")

    # Second expansion
    t3_term = t2_term * 2
    c3_val = t2_term * c1_val + c2_val
    print("\nAfter substituting T(n-2) = 2*T(n-3) + c:")
    print(f"  T(n) = {t2_term} * (2*T(n-3) + c) + {c2_val}*c")
    print(f"  T(n) = {t3_term} * T(n - 3) + {c3_val} * c")
    
    print("\nWe see a pattern emerging: T(n) = 2^k * T(n-k) + (2^k - 1) * c")
    print("The recursion stops at the base case T(1), which is a constant. This occurs when n - k = 1, so k = n - 1.")
    
    print("\nSubstituting k = n - 1 gives the final equation for T(n):")
    print("  T(n) = 2^(n - 1) * T(1) + (2^(n - 1) - 1) * c")
    
    print("\nThe dominant term in this expression is 2^(n-1), as T(1) and c are constants.")
    print("Therefore, the time complexity is O(2^(n-1)), which simplifies to O(2^n).")

if __name__ == '__main__':
    main()