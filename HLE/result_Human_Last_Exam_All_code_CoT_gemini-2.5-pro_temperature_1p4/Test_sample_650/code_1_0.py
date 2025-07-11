def explain_complexity():
    """
    Explains the computational time complexity of the IMAG(n, c) algorithm.
    """
    print("Computational Time Complexity Analysis of IMAG(n, c):")
    print("-" * 55)
    
    print("1. The algorithm's execution time is primarily determined by the 'while q > 0' loop.")
    
    print("\n2. Let's analyze the number of loop iterations. The loop continues as long as q > 0.")
    print("   - Initially, q is set to floor(n / c).")
    print("   - In each loop iteration, the effective value of q is divided by c.")
    print("   - The sequence of values for q is approximately n/c, n/(c^2), n/(c^3), ..., until the value is less than 1 (making q=0).")
    
    print("\n3. The number of times you can divide 'n' by 'c' until the result is less than 1 is the definition of log_c(n).")
    print("   - Therefore, the loop runs approximately log_c(n) times.")

    print("\n4. The operations inside the loop (incrementing i, assignment, division, multiplication, subtraction) are all constant-time operations.")
    
    print("\n5. The total time complexity is the (number of iterations) * (time per iteration).")
    print(f"   - Complexity = O(log_c(n)) * O(1) = O(log_c(n))")

    print("\n6. In Big O notation, the base of the logarithm is usually omitted because it only differs by a constant factor (log_c(n) = log(n) / log(c)).")
    print("   - The final complexity is simplified to O(log n).")

    print("\nFinal Conclusion:")
    print("The time complexity of the algorithm is O(log_c(n)), which is commonly expressed as O(log n).")

explain_complexity()