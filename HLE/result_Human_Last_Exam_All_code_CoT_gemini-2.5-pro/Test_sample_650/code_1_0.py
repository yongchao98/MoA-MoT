import math

def analyze_imag_complexity():
    """
    This function implements the IMAG algorithm to analyze and explain
    its computational time complexity. It also constructs the final
    equation representing the base conversion.
    """

    # 1. Define the algorithm as a Python function
    def IMAG(n, c):
        """Converts an integer n to its base c representation."""
        if not (isinstance(n, int) and isinstance(c, int) and n >= 0 and c >= 2):
            raise ValueError("Inputs must be integers where n >= 0 and c >= 2.")
        
        if n == 0:
            return [0]

        a = []
        x = n
        # Step 1 of the algorithm
        q = x // c
        a.append(x - q * c)
        
        # Step 2 of the algorithm (the while loop)
        while q > 0:
            x = q
            q = x // c
            a.append(x - q * c)
            
        # The algorithm produces digits in reverse order (a_0, a_1, ...), 
        # so we reverse the list for the standard representation (a_i, a_{i-1}, ...).
        return a[::-1]

    # --- Complexity Analysis ---
    
    # We will use example values to make the explanation concrete.
    n_val = 987
    c_val = 8

    # Run the algorithm to get the base-c representation digits (the a_i values).
    representation = IMAG(n_val, c_val)
    num_digits = len(representation)

    # Print the step-by-step complexity analysis.
    print("Computational Time Complexity Analysis of IMAG(n, c):")
    print("-" * 60)
    print("1. The algorithm's runtime is dominated by its 'while' loop.")
    print("2. The loop runs as long as the quotient 'q' is greater than 0.")
    print("3. In each iteration, the value of 'q' is reduced by a factor of 'c'.")
    print(f"4. The number of loop iterations is therefore O(log_c(n)). This is proportional to the number of digits in the output.")
    print(f"   - For n={n_val}, c={c_val}, the output has {num_digits} digits.")
    print(f"   - The theoretical number of iterations is floor(log_{c_val}({n_val})) = {math.floor(math.log(n_val, c_val))}.")
    print("5. The operations inside the loop (division, multiplication, etc.) take constant time, O(1).")
    print("6. Total Complexity = (Iterations) * (Work per iteration) = O(log_c(n)) * O(1).")
    print("\nFinal Complexity: O(log_c(n))")
    print("-" * 60)

    # The "final equation" is the mathematical definition of a base representation,
    # which the algorithm calculates. We print it here.
    print("The algorithm finds the coefficients 'a_i' for the equation:")
    print("n = a_i*c^i + a_{i-1}*c^{i-1} + ... + a_0*c^0\n")
    
    equation_parts = []
    # Iterate through the digits from most significant (a_i) to least significant (a_0)
    for i, digit in enumerate(representation):
        power = num_digits - 1 - i
        # Add each number and term to the equation string
        equation_parts.append(f"{digit} * {c_val}^{power}")
        
    equation_str = " + ".join(equation_parts)
    print(f"For the given inputs, the final equation with each number is:")
    print(f"{n_val} = {equation_str}")

# Execute the analysis function.
analyze_imag_complexity()
<<<O(log_c n)>>>