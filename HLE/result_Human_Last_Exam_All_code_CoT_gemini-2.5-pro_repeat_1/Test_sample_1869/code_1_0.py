import sys

# In Python 3.11+, the default limit for int to string conversion
# might be too low for the very large numbers in this calculation.
# We can increase it if needed, though it's often not necessary
# for the calculation itself, just for printing.
# sys.set_int_max_str_digits(0) # 0 means no limit

def find_special_prime():
    """
    Finds the smallest prime p from a list of choices such that
    the ring Z[p√6] is not the ring of integers of the number field Q(p√6).

    This condition is equivalent to finding a prime p such that 6^(p-1) ≡ 1 (mod p^2).
    """
    
    print("The condition for Z[p√6] to not be the ring of integers of Q(p√6) for a prime p > 3 is:")
    print("6^(p-1) ≡ 1 (mod p^2)\n")

    answer_choices = {
        17: 'A',
        383: 'B',
        1093: 'C',
        66161: 'D',
        534851: 'E'
    }
    
    primes = sorted(answer_choices.keys())
    
    found = False
    for p in primes:
        print(f"--- Checking prime p = {p} (Choice {answer_choices[p]}) ---")
        
        p_squared = p * p
        
        # Use pow(base, exp, mod) for efficient modular exponentiation
        result = pow(6, p - 1, p_squared)
        
        print(f"Calculating 6^({p-1}) mod {p}^2...")
        print(f"The result is: {result}")
        
        if result == 1:
            print(f"\nCondition met! For p = {p}, we have 6^({p-1}) ≡ 1 (mod {p}^2).")
            print(f"This is the smallest prime in the list that satisfies the condition.")
            
            print("\nTo output the full equation 6^(p-1) = k * p^2 + 1, we need to find k.")
            
            # Use Python's arbitrary-precision integers to calculate k
            # k = (6^(p-1) - 1) / p^2
            # This might take a few seconds as the numbers are very large.
            print("Calculating k... (this may take a moment)")
            try:
                # Calculate the large power first
                large_power = pow(6, p - 1)
                # Then calculate k
                k = (large_power - 1) // p_squared
                
                print(f"The final equation is:\n6^({p-1}) = {k} * {p}^2 + 1")

            except Exception as e:
                print(f"Could not compute k due to an error: {e}")
                print(f"However, the prime is correctly identified as {p}.")

            found = True
            break
        else:
            print(f"Condition not met for p = {p}, since the result is not 1.\n")
            
    if not found:
        print("No prime in the list satisfies the condition.")

if __name__ == '__main__':
    find_special_prime()