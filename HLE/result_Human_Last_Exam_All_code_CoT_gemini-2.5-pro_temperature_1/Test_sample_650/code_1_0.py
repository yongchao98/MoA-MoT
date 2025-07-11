import math

def analyze_and_run_IMAG(n, c):
    """
    This function implements the IMAG algorithm to convert n to its base c representation.
    It also prints the final representation, a verification equation, and data for 
    complexity analysis.
    
    Args:
        n (int): The integer to convert (n >= 0).
        c (int): The base to convert to (c >= 2).
    """
    # --- Start of the IMAG Algorithm Implementation ---
    
    if n < 0 or c < 2:
        print("Input error: n must be >= 0 and c must be >= 2.")
        return

    print(f"Running IMAG(n={n}, c={c})")
    
    if n == 0:
        print("The base", c, "representation of", n, "is: (0)")
        print("Final Equation: 0 = 0 *", c, "^ 0")
        print("Number of loop iterations: 0")
        return

    # Step 1: Initialization
    i = 0
    x = n
    q = x // c
    a0 = x - q * c
    
    # List 'a' stores digits in reverse order: a0, a1, ...
    a = [a0]
    
    loop_iterations = 0

    # Step 2: Loop
    while q > 0:
        loop_iterations += 1
        # Step 2.1
        i += 1
        x = q
        q = x // c
        ai = x - q * c
        a.append(ai)
        
    # Step 3: Final result
    # Reverse the list to get the correct representation (ai, ..., a1, a0)
    representation = a[::-1]
    
    # --- End of Algorithm Implementation ---


    # --- Output and Analysis ---
    
    # Print the representation as a tuple
    print("The base", c, "representation of", n, "is:", tuple(representation))

    # Print the final equation showing how each digit contributes to the original number
    equation_parts = []
    check_sum = 0
    power = len(representation) - 1
    for digit in representation:
        equation_parts.append(f"{digit} * {c}^{power}")
        check_sum += digit * (c ** power)
        power -= 1
    
    print("Final Equation:", n, "=", " + ".join(equation_parts))
    print(f"Verification: The sum evaluates to {check_sum}.")
    print("-" * 20)
    
    # Print the complexity analysis based on the run
    print("Complexity Analysis:")
    print(f"The 'while' loop ran {loop_iterations} times.")
    
    # The theoretical number of iterations is floor(log_c(n))
    theoretical_iterations = math.floor(math.log(n, c))
    print(f"Theoretical number of iterations (floor(log_c(n))) is: {theoretical_iterations}")
    print("This shows that the number of operations is proportional to log(n).")

# --- Main execution ---
# You can change these values to test the algorithm
input_n = 2024
input_c = 7

analyze_and_run_IMAG(input_n, input_c)