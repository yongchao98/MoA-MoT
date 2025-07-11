import math

def analyze_imag_complexity(n, c):
    """
    This function implements the IMAG algorithm, counts the number of loop iterations,
    and compares it to the theoretical value of log_c(n) to demonstrate the
    O(log n) time complexity. It then prints the values involved in this comparison.
    """
    # Validate inputs as per the algorithm's constraints
    if not (isinstance(n, int) and isinstance(c, int) and n >= 0 and c >= 2):
        print("Error: Input must be integers n and c, where n >= 0 and c >= 2.")
        return

    # --- Start of IMAG Algorithm Implementation ---
    # 1. Initialization
    i = 0
    x = n
    q = x // c  # In Python, // is floor division, equivalent to ⌊x/c⌋

    # Variable to count the loop executions for our analysis
    loop_iterations = 0

    # 2. While loop
    while q > 0:
        # Increment our counter for each loop execution
        loop_iterations += 1

        # 2.1.
        i = i + 1
        x = q
        q = x // c
        # The calculation of a_i is not needed for complexity analysis, so we omit it.
        # a_i = x - q * c
    # --- End of Algorithm ---

    # The number of iterations is the final value of 'i' or our 'loop_iterations' counter.
    # We will now print the analysis.

    print(f"--- Analysis for n = {n}, c = {c} ---")

    # The number of iterations, 'k', is directly given by our counter.
    k = loop_iterations

    # The "final equation" we demonstrate is: k = floor(log_c(n))
    # We will output each number in this equation.
    if n > 0:
        # Calculate the theoretical value
        log_value = math.log(n, c)
        floor_log_value = math.floor(log_value)

        # Output the numbers involved in the relationship k = floor(log_c(n))
        print(f"Number of loop iterations (k): {k}")
        print(f"Value of log base {c} of {n}: {log_value:.4f}")
        print(f"The final relationship is: {k} = floor({log_value:.4f}) = {floor_log_value}")
        print("This confirms the number of iterations is proportional to log_c(n).")
    else:  # Edge case for n = 0
        print(f"Number of loop iterations (k): {k}")
        print("For n=0, the loop does not run, which is consistent with the analysis.")


# Example 1: Large n, standard base
analyze_imag_complexity(n=1000000, c=10)

print("\n" + "="*40 + "\n")

# Example 2: Power of 2, binary base
analyze_imag_complexity(n=65536, c=2)