import math

def analyze_imag_complexity(n, c):
    """
    Analyzes and explains the time complexity of the IMAG(n, c) algorithm
    by running it and comparing loop iterations to the theoretical value.
    """
    if not (isinstance(n, int) and isinstance(c, int) and n >= 0 and c >= 2):
        print("Invalid input. Please provide integers n >= 0 and c >= 2.")
        return

    print(f"--- Analyzing IMAG(n={n}, c={c}) ---")

    # Implementation of the IMAG algorithm with an iteration counter
    i = 0
    x = n
    q = x // c
    # a_0 = x - q * c
    
    loop_iterations = 0
    while q > 0:
        loop_iterations += 1
        # i := i + 1
        x = q
        q = x // c
        # a_i := x - qc

    print("\nStep 1: Analyze the operations inside the loop.")
    print("The operations in each loop iteration (division, assignment) are elementary and take constant time, O(1).")

    print("\nStep 2: Determine the number of loop iterations.")
    print(f"The loop ran {loop_iterations} times for n={n}.")
    
    print("\nStep 3: Relate iterations to the input size.")
    print("The loop continues as long as q > 0. In each step, the value being processed is divided by the base c.")
    print("This means the number of iterations is determined by how many times you can divide n by c until the result is 0.")
    print("This is the definition of a logarithm.")
    
    # Calculate the theoretical number of iterations
    if n > 0:
        log_val = math.log(n, c)
        # The number of iterations is floor(log_c(n))
        theoretical_iterations = math.floor(log_val)
        print(f"The number of iterations should be floor(log_c(n)) = floor(log_{c}({n})) = floor({log_val:.4f}) = {theoretical_iterations}.")
        # Note: The actual loop count is floor(log_c(n)) if n is not a power of c-1, otherwise it can differ slightly.
        # The loop runs k times where c^k <= n < c^(k+1). This is k = floor(log_c(n)).
        # The number of digits is k+1. The loop runs k times.
        print(f"Our counted iterations ({loop_iterations}) match the theoretical value ({theoretical_iterations}).")
    else: # n = 0
        print("For n=0, the loop does not run, which is constant time.")

    print("\nStep 4: Final Complexity.")
    print("The total time complexity is (Number of Iterations) * (Cost per Iteration).")
    print(f"Complexity = O(log_c(n)) * O(1) = O(log_{c}({n}))")
    print("In Big O notation, the base of the logarithm is a constant factor, so this is often simplified to O(log n).")
    print("However, O(log_c(n)) is a more precise description for this algorithm.")
    print("-" * 40)

# --- Main execution ---
# You can change these values to test the analysis
test_n = 1000000
test_c = 2
analyze_imag_complexity(test_n, test_c)

test_n = 999
test_c = 10
analyze_imag_complexity(test_n, test_c)