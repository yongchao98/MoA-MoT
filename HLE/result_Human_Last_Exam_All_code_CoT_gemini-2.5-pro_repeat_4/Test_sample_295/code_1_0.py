import sys

def compute_euler_char_mod_k(k):
    """
    This function computes the reduced Euler characteristic of Delta_k modulo k.

    Args:
        k: A prime number k >= 3.
    """
    # First, we check if k is a prime number greater than or equal to 3.
    if not isinstance(k, int) or k < 3:
        print(f"Error: Input k must be an integer >= 3. Received: {k}", file=sys.stderr)
        return
    
    is_prime = True
    if k > 3:
        for i in range(2, int(k**0.5) + 1):
            if k % i == 0:
                is_prime = False
                break
    
    if not is_prime:
        print(f"Error: Input k must be a prime number. {k} is not prime.", file=sys.stderr)
        return

    # Based on the derivation, the formula for the result is (k-1)/2.
    # The problem asks to output each number in the final equation.
    
    # The numbers in the equation
    k_val = k
    one = 1
    two = 2
    
    # Perform the calculation
    numerator = k_val - one
    result = numerator // two
    
    # The result modulo k is the result itself, as 0 <= (k-1)/2 < k.
    final_result_mod_k = result % k
    
    print(f"For the prime k = {k_val}:")
    print(f"The formula for the reduced Euler characteristic modulo k is (k - 1) / 2.")
    print(f"Step 1: Calculate the numerator (k - 1)")
    print(f"   {k_val} - {one} = {numerator}")
    print(f"Step 2: Calculate the division")
    print(f"   {numerator} / {two} = {result}")
    print(f"The final result is {result}.")
    print(f"Checking the result modulo k: {result} mod {k_val} = {final_result_mod_k}.")


# We will use an example prime k=19 for demonstration.
example_k = 19
compute_euler_char_mod_k(example_k)
