import sys
from mpmath import mp
import sympy

def solve_irrational_prime_puzzle():
    """
    Solves the puzzle by searching through common irrational numbers
    for a specific prime sequence property.
    """
    # Set the precision for the mpmath library to 500 digits.
    mp.dps = 500
    # The maximum length of the number to check for primality.
    # Based on known sequences, the 6th prime is likely found within 400 digits.
    MAX_DIGIT_LENGTH = 450

    # Define the candidate irrational numbers and their symbols.
    # The problem specifies a symbol of no more than 3 characters (e.g., 'e', 'π', '√2').
    candidates = {
        "π": mp.pi,
        "e": mp.e,
        "φ": mp.phi,
        "√2": mp.sqrt(2),
    }

    # Iterate through each candidate irrational number.
    for symbol, num_obj in candidates.items():
        print(f"--- Checking irrational number: {symbol} ---")
        
        # Get the digits of the number as a string, removing the decimal point.
        digit_string = str(num_obj).replace('.', '')
        
        primes_found = []
        
        # Check prefixes of increasing length for primality.
        for length in range(1, MAX_DIGIT_LENGTH + 1):
            if length > len(digit_string):
                break

            # Create the number from the prefix.
            prefix_str = digit_string[:length]
            try:
                candidate_num = int(prefix_str)
            except (ValueError, TypeError):
                continue
                
            # Check if the number is prime.
            if sympy.isprime(candidate_num):
                primes_found.append(candidate_num)
                
                # If we have found the 6th prime, check its last digits.
                if len(primes_found) == 6:
                    sixth_prime = primes_found[-1]
                    # Check if the last 6 digits match the target "521023".
                    if str(sixth_prime).endswith("521023"):
                        print(f"\n>>> Match Found for irrational number: {symbol} <<<")
                        print("The first 6 primes generated from its digits are:")
                        for i, p in enumerate(primes_found, 1):
                            p_str = str(p)
                            if len(p_str) > 20:
                                p_str = f"{p_str[:10]}...{p_str[-10:]}"
                            print(f"  {i}: {p_str} (length {len(str(p))})")
                        
                        print("\nThe 6th prime number in the sequence is:")
                        # The problem asks to output the numbers in the final equation.
                        # We print the full 6th prime to satisfy this.
                        print(str(sixth_prime))
                        
                        # Return the symbol of the found irrational number.
                        return symbol
                    else:
                        # If the 6th prime doesn't match, this isn't our number.
                        # We can stop checking this irrational and move to the next.
                        print("Found 6 primes, but the last one does not match. Moving to next candidate.")
                        break # Exit the inner loop and go to the next irrational number
                        
        print(f"Finished checking {symbol}. No match found within {MAX_DIGIT_LENGTH} digits.\n")

    return "No solution found among the candidates."

if __name__ == '__main__':
    solution = solve_irrational_prime_puzzle()
    print(f"\nFinal Answer: The irrational number is identified by the symbol {solution}.")
    sys.stdout.flush()
    # The final answer format required by the user prompt.
    print(f"<<<{solution}>>>")
