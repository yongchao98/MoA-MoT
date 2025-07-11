import decimal
import random

def is_prime(n):
    """
    Checks if a number n is prime using the Miller-Rabin primality test.
    This method is efficient for large numbers.
    """
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    
    # Write n as 2^r * d + 1
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1

    # Perform k tests for accuracy. k=5 is sufficient for our purpose.
    for _ in range(5):
        a = random.randint(2, n - 2)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
            
    return True

def find_special_prime():
    """
    Calculates digits of 'e', finds the 6th prime in its concatenated
    digit sequence, and checks if it matches the puzzle's condition.
    """
    # Set precision for decimal calculations. 250 digits should be more than enough.
    decimal.getcontext().prec = 250
    
    # Calculate Euler's number 'e' and create the digit string
    e_digits = str(decimal.Decimal(1).exp()).replace('.', '')
    
    target_ending = "521023"
    primes_found = []
    
    # Iterate through prefixes of the digit string
    for i in range(1, len(e_digits) + 1):
        num_str = e_digits[:i]
        num = int(num_str)
        
        # A small optimization: quickly discard numbers ending in 0, 2, 4, 5, 6, 8
        if i > 1 and num % 10 in (0, 2, 4, 5, 6, 8):
            continue

        if is_prime(num):
            primes_found.append(num)
            if len(primes_found) == 6:
                sixth_prime = num
                # Check if we found the solution
                if str(sixth_prime).endswith(target_ending):
                    print("Solution Found!")
                    print("The irrational number is Euler's number, expressed with the symbol 'e'.")
                    print("\nThe 6th prime number generated from its digits is:")
                    print(sixth_prime)
                    print(f"\nThis number ends in {str(sixth_prime)[-6:]}, which matches the required {target_ending}.")
                break

# Execute the function to solve the puzzle
find_special_prime()
