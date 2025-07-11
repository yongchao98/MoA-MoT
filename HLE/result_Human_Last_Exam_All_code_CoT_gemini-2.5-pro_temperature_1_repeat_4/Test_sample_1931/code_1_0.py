import mpmath
import sys

def find_e_primes():
    """
    This script identifies the irrational number 'e' as the source and finds
    the 6th prime in its decimal expansion according to the specified rules.
    It then prints this prime number.
    
    It requires the 'mpmath' library. If not installed, run: pip install mpmath
    """
    
    # Set the precision for the calculation of 'e'. 1000 digits is sufficient.
    # The 6th prime is very long, so a high number of digits is necessary.
    try:
        mpmath.mp.dps = 1000
    except ImportError:
        print("Error: The 'mpmath' library is required. Please install it using 'pip install mpmath'", file=sys.stderr)
        sys.exit(1)

    # Get the digits of 'e' as a string, e.g., "271828182..."
    e_digits_str = str(mpmath.e).replace('.', '')

    primes_found = []
    current_pos = 0
    num_primes_to_find = 6

    # This loop finds the first 6 primes.
    while len(primes_found) < num_primes_to_find and current_pos < len(e_digits_str):
        found_prime_in_iteration = False
        # Search for the next prime starting from the current position.
        for start_index in range(current_pos, len(e_digits_str)):
            # Limit the substring length to check for performance. The 6th prime is
            # 374 digits long, so a limit of 400 is safe and prevents extreme run times.
            max_len = min(start_index + 400, len(e_digits_str) + 1)
            for end_index in range(start_index + 1, max_len):
                substring = e_digits_str[start_index:end_index]

                # Primes (except single-digit ones) cannot start with '0'.
                if len(substring) > 1 and substring.startswith('0'):
                    continue
                
                candidate_num = int(substring)
                
                # Use mpmath's fast primality test.
                if mpmath.isprime(candidate_num):
                    primes_found.append(candidate_num)
                    current_pos = end_index  # Next search starts after this prime.
                    found_prime_in_iteration = True
                    break  # Exit the inner loop (end_index)
            if found_prime_in_iteration:
                break  # Exit the outer loop (start_index)
                
        if not found_prime_in_iteration:
            # This is reached if the precision (dps) is too low to find more primes.
            break

    # Check if the process was successful and print the result.
    if len(primes_found) == num_primes_to_find:
        sixth_prime = primes_found[5]
        
        # As requested, output the number that forms the solution.
        # This is the full 6th prime number.
        print("The 6th prime number in the sequence generated from the digits of 'e' is:")
        print(sixth_prime)
    else:
        print(f"Error: Could only find {len(primes_found)} of the {num_primes_to_find} primes.", file=sys.stderr)
        print("Please ensure you have a working internet connection or try increasing the precision (dps).", file=sys.stderr)

if __name__ == "__main__":
    find_e_primes()