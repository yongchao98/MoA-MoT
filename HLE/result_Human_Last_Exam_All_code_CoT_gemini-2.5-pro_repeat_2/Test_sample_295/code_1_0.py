def compute_reduced_euler_characteristic_mod_k(k):
    """
    Computes the value of the reduced Euler characteristic of the complex Delta_k,
    modulo k, for a given prime k >= 3.

    The theoretical derivation shows that the result is (k-1)/2.
    This function implements this formula and prints the steps.
    
    Args:
        k (int): A prime number such that k >= 3.
    """
    if k < 3 or (k > 2 and k % 2 == 0):
        print("Error: k must be a prime number greater than or equal to 3.")
        return

    # The formula derived from the group-theoretic argument is (k-1)/2.
    # We will show the numbers in the equation.
    
    print(f"The problem is to compute the reduced Euler characteristic of Delta_k mod k.")
    print(f"For k = {k}, we use the derived formula.")
    
    # The final formula for the result mod k is (k-1)/2
    final_formula = "(k - 1) / 2"
    print(f"The formula for the result is: {final_formula}")
    
    # Numerator of the formula
    numerator = k - 1
    print(f"Numerator: k - 1 = {k} - 1 = {numerator}")
    
    # Denominator of the formula
    denominator = 2
    print(f"Denominator: {denominator}")

    # Since k is an odd prime, k-1 is even, so the division results in an integer.
    result = numerator // denominator
    
    # The result (k-1)/2 is an integer smaller than k, so the modulo operation is trivial.
    print(f"Result: {numerator} / {denominator} = {result}")
    print(f"Therefore, the value of the reduced Euler characteristic of Delta_{k} mod {k} is {result}.")
    
    # For the final answer format
    global final_answer
    final_answer = result


# --- Main execution ---
# Let's choose a prime k >= 3 as an example. The user can change this value.
# For instance, let's use k = 17.
k_value = 17 
final_answer = None
compute_reduced_euler_characteristic_mod_k(k_value)
# The final result for k=17 is printed and stored in final_answer.
# The following will not be printed in a standard run but is here to generate the final answer block.
# print(f"<<<{final_answer}>>>")