import mpmath
import random

def is_prime(n):
    """
    Checks if a number n is prime using the Miller-Rabin primality test.
    This test is probabilistic but highly accurate for a sufficient number of rounds.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False

    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # Run the test 5 times for high confidence
    for _ in range(5):
        a = random.randrange(2, n - 1)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def find_irrational_number():
    """
    Finds the irrational number that matches the puzzle's criteria.
    """
    # Set mpmath precision. We need ~350 digits for the 6th prime. 400 is safe.
    mpmath.mp.dps = 400
    
    # Generate the digits of Euler's number 'e' and remove the decimal point.
    digit_string = str(mpmath.e).replace('.', '')
    
    primes_found = []
    print("Searching for primes in the initial digits of e (2.71828...)...\n")
    
    # Iterate through prefixes of the digit string
    for length in range(1, len(digit_string) + 1):
        num_str = digit_string[0:length]
        
        # Convert string to integer for primality testing
        n = int(num_str)
        
        if is_prime(n):
            primes_found.append(n)
            print(f"Found prime #{len(primes_found)} (length {len(num_str)}): {n}")
            
            # Check if we have found the 6th prime
            if len(primes_found) == 6:
                print("\nThis is the 6th prime. Let's check its last 6 digits.")
                sixth_prime_str = str(primes_found[5])
                
                if sixth_prime_str.endswith("521023"):
                    print(f"\nSUCCESS: The last 6 digits are indeed '521023'.")
                    print("The irrational number is Euler's number, e.")
                    # Per the instruction to "output each number in the final equation!"
                    # we print the full 6th prime.
                    print("\nThe final prime number is:")
                    print(sixth_prime_str)
                else:
                    print(f"The last 6 digits ({sixth_prime_str[-6:]}) do not match.")
                
                # Stop the search
                return
                
    print("\nCould not find the 6th prime within the calculated precision.")

if __name__ == '__main__':
    find_irrational_number()