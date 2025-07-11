import math

def analyze_imag_complexity(n, c):
    """
    This function analyzes the IMAG algorithm by counting loop iterations
    and printing the final complexity equation.
    """
    if not (isinstance(n, int) and isinstance(c, int) and n >= 0 and c >= 2):
        print("Invalid input: Please provide integers n >= 0 and c >= 2.")
        return

    print(f"--- Analyzing IMAG algorithm for n={n}, c={c} ---")

    # Initialization (O(1) time)
    i = 0
    x = n
    q = x // c
    
    loop_iterations = 0
    
    # The 'while' loop is the dominant part of the algorithm.
    while q > 0:
        loop_iterations += 1
        # Operations inside the loop are considered O(1).
        x = q
        q = x // c

    print(f"The 'while' loop executed {loop_iterations} times.")

    # The number of iterations is proportional to log_c(n).
    # We can verify this calculation.
    if n >= c:
        theoretical_log = math.log(n, c)
        theoretical_iterations = math.floor(theoretical_log)
        print(f"The theoretical number of iterations is floor(log_c(n)) = floor(log_{c}({n})) = floor({theoretical_log:.4f}) = {theoretical_iterations}.")
    else:
        print("The loop does not run as n < c (0 iterations).")

    # Final Complexity Equation
    print("\nThe computational time complexity is the number of iterations multiplied by the constant work per iteration.")
    print("This gives the final equation: T(n, c) = O(log_c(n))")
    print("Since c is a constant, this is equivalent to T(n) = O(log n).")
    print("-" * 50)


# --- Main execution with example values ---
analyze_imag_complexity(1000, 2)
analyze_imag_complexity(1000000, 10)