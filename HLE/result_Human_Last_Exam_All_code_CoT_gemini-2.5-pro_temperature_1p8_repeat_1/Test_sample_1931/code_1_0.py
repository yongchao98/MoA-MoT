import sympy
import sys

def solve_puzzle():
    """
    Finds the irrational number that satisfies the puzzle's conditions.
    This script requires the 'sympy' library. You can install it by running:
    pip install sympy
    """
    # The last 6 digits we are looking for in the 6th prime.
    target_last_digits = 521023
    
    # We need to find the 6th prime.
    target_prime_index = 6

    # Set the precision for the calculation of irrational numbers' digits.
    # The 6th prime is quite long, so we need a high precision. 400 digits is sufficient.
    precision = 400

    # Define the candidate irrational numbers.
    # The name, the sympy object, and the number of symbols for its expression.
    candidates = [
        ("e", sympy.E, 1),
        ("pi", sympy.pi, 1),
        ("the Golden Ratio (phi)", sympy.GoldenRatio, 1),
        ("the square root of 2", sympy.sqrt(2), 2),
    ]

    print("Starting the search for the irrational number...")
    print("----------------------------------------------")

    for name, number_obj, _ in candidates:
        print(f"\n>>> Testing {name} <<<")
        
        # Get the digits of the irrational number as a string, removing the decimal point.
        digit_str = str(sympy.N(number_obj, precision)).replace('.', '')

        primes_found = []
        current_num_str = ''

        for digit in digit_str:
            current_num_str += digit

            # Basic primality test optimizations:
            # - Primes > 5 must end in 1, 3, 7, or 9.
            # - '1' is not a prime number.
            if len(current_num_str) > 1 and current_num_str[-1] in '024568':
                continue
            if current_num_str == '1':
                continue

            num = int(current_num_str)

            # Use sympy's robust primality test.
            if sympy.isprime(num):
                primes_found.append(num)

            # Check if we have found the 6th prime.
            if len(primes_found) == target_prime_index:
                sixth_prime = primes_found[-1]
                
                # Check if the last 6 digits match the target.
                if sixth_prime % 10**6 == target_last_digits:
                    print(f"\nSUCCESS! The irrational number is {name}.")
                    print("\nThe problem is solved.")
                    print("\nHere are the details:")
                    print(f"\nThe 6th prime number generated from the digits of {name} is:")
                    print(f"\n{sixth_prime}")
                    
                    print(f"\nTo confirm, we check the final condition:")
                    print(f"The number {sixth_prime} modulo 1,000,000 is {sixth_prime % 10**6}.")
                    
                    print(f"\nFinal Equation:")
                    print(f"{sixth_prime} % 10**6 = {target_last_digits}")
                    
                    print(f"\nThis shows that {target_last_digits} are indeed the last 6 digits.")
                    return name
                
                # If the digits don't match, this is not our irrational number.
                # We can stop checking it and move to the next candidate.
                else:
                    print(f"Found the 6th prime, but its last digits do not match. Moving to the next candidate.")
                    break # Breaks from the inner digit-loop to the next candidate
        
        # This message prints if we exhausted the precision without finding 6 primes.
        if len(primes_found) < target_prime_index:
            print(f"Could not find {target_prime_index} primes within the first {precision} digits.")

    print("\nSearch complete. No matching irrational number found among the candidates.")
    return None

if __name__ == '__main__':
    # Ensure sympy is installed before running.
    try:
        import sympy
    except ImportError:
        print("Error: The 'sympy' library is required to run this script.")
        print("Please install it using: pip install sympy")
        sys.exit(1)
        
    solve_puzzle()
