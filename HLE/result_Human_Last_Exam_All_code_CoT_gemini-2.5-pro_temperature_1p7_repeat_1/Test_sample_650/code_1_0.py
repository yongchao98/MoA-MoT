import math

def print_complexity_analysis():
    """
    This function prints a step-by-step analysis of the IMAG(n, c) algorithm's
    computational time complexity.
    """
    print("Computational Time Complexity Analysis of IMAG(n, c):")
    print("----------------------------------------------------------")
    print("1. The algorithm consists of an initialization step and a 'while' loop.")
    print("   The initialization takes a constant amount of time, O(1).")
    print("   The algorithm's overall complexity is determined by the 'while' loop.")
    print("\n2. The 'while' loop continues as long as q > 0. Let's see how 'q' changes:")
    print("   - Before the loop: q = floor(n / c)")
    print("   - After 1st iteration: q = floor(n / c^2)")
    print("   - After k-th iteration: q = floor(n / c^(k+1))")
    print("\n3. The loop stops when q = 0, which happens when n / c^(k+1) < 1, or c^(k+1) > n.")
    print("\n4. To find the number of iterations k, we can solve for k:")
    print("   log_c(c^(k+1)) > log_c(n)")
    print("   k + 1 > log_c(n)")
    print("   k > log_c(n) - 1")
    print("   This shows that the number of iterations k is proportional to log_c(n).")
    print("\n5. The operations inside the loop (addition, division, assignment) take constant time, O(1).")
    print("\n6. Therefore, the total time complexity T(n) is the number of iterations times the constant work per iteration.")
    print("\n   T(n) = (number of iterations) * (work per iteration)")
    print("   T(n) = O(log_c(n)) * O(1)")
    print("\nFinal Equation for Complexity:")
    # The base of the logarithm, c, is a number (integer >= 2) from the algorithm's input.
    # The final equation includes this number 'c' as the base.
    print("   T(n) = O(log_c n)")

print_complexity_analysis()