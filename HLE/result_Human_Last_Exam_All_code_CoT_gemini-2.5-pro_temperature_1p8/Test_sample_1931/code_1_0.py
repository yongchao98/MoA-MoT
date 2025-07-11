import mpmath
from sympy import isprime

def find_special_prime():
    """
    This program solves the puzzle by identifying an irrational number and
    the specific prime in its digit sequence that matches the given criteria.

    The puzzle requires finding an irrational number, expressible with no
    more than 3 symbols, where the 6th prime number generated from its
    concatenated digits ends in 521023.

    The chosen candidate is 1/e, which satisfies the symbol constraint.
    """
    
    # Set the required precision. The prime number we are looking for has 36 digits,
    # so a precision of 50 is sufficient.
    mpmath.mp.dps = 50
    
    # Define the target properties
    target_prime_index = 6
    target_ending = 521023
    
    # Generate the digit string for our candidate, 1/e.
    # 1/e = 0.36787944...
    # We take the digits after the decimal point: "36787944..."
    irrational_name = "1/e"
    try:
        val = 1 / mpmath.e
        digit_string = str(val).replace('.', '').lstrip('0')
    except ImportError:
        print("Please install mpmath and sympy: pip install mpmath sympy")
        return
    except Exception as e:
        print(f"An error occurred: {e}")
        return

    prime_count = 0
    final_prime = None
    
    print(f"Searching for the {target_prime_index}th prime in the digit sequence of {irrational_name}...")

    # Iterate through prefixes of the digit string to form numbers
    for k in range(1, len(digit_string) + 1):
        num_str = digit_string[:k]
        
        # An empty string would cause an error, so we skip.
        if not num_str:
            continue
        
        num = int(num_str)
        
        # Primes must be greater than 1.
        if num < 2:
            continue
            
        # Test for primality
        if isprime(num):
            prime_count += 1
            
            # Check if this is the 6th prime number
            if prime_count == target_prime_index:
                final_prime = num
                # Found what we are looking for, so we can exit the loop.
                break

    # After the search, present the results.
    if final_prime is not None:
        print("\n--- Solution Found ---")
        print(f"The irrational number is: {irrational_name}")
        
        # The prompt asks to output each number in the "final equation"
        print(f"The 6th prime number found in the sequence is:")
        print(final_prime)
        
        last_six_digits = final_prime % 1000000
        print(f"\nVerifying the last 6 digits of this prime:")
        print(f"{final_prime} % 1,000,000 = {last_six_digits}")

        if last_six_digits == target_ending:
            print(f"\nThe result {last_six_digits} matches the target {target_ending}.")
        else:
            print(f"The result {last_six_digits} does NOT match the target {target_ending}.")
    else:
        print("Failed to find the 6th prime within the calculated precision.")


if __name__ == "__main__":
    find_special_prime()