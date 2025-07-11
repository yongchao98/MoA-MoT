import math

def analyze_imag_complexity(n, c):
    """
    Analyzes the computational time complexity of the IMAG(n, c) algorithm.
    It counts the loop iterations required to convert n to base c and prints
    the complexity analysis.
    """
    if not (isinstance(n, int) and isinstance(c, int) and n >= 0 and c >= 2):
        print("Error: Input must be integers n, c, where n >= 0 and c >= 2.")
        return

    # Let's trace the algorithm to count the loop iterations.
    x = n
    q = x // c
    
    loop_iterations = 0
    
    # The 'while' loop is the core of the algorithm's complexity.
    while q > 0:
        loop_iterations += 1
        x = q
        q = x // c
        
    print(f"Analysis for IMAG(n={n}, c={c}):")
    print(f"The 'while' loop executed {loop_iterations} times.")
    
    # The number of iterations is related to log_c(n).
    # The loop runs floor(log_c(n)) times for n >= c.
    if n >= c:
        log_c_n = math.log(n, c)
        theoretical_iterations = math.floor(log_c_n)
        print(f"The value of log_{c}({n}) is approximately {log_c_n:.4f}.")
        print(f"The theoretical number of iterations, floor(log_c(n)), is {theoretical_iterations}.")
    else: # Case where 0 <= n < c
        print("Since n < c, the initial quotient 'q' is 0, so the loop does not run.")

    print("\n--- Complexity Conclusion ---")
    print("The number of loop iterations is proportional to log_c(n).")
    print("Since the operations inside the loop are constant time (O(1)), the total time complexity is determined by the number of iterations.")
    print("The final complexity equation is: T(n) = O(log n)")
    
    # Printing each component of the final equation
    print("\nFinal Equation Components:")
    print("O", "(", "log", "n", ")")

# --- Example Execution ---
# You can change these values to test different scenarios.
n_value = 1000
c_value = 2
analyze_imag_complexity(n_value, c_value)