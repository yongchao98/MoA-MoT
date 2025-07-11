import decimal
import math

def is_prime(n):
    """
    Checks if a number n is prime using an optimized trial division method.
    """
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_irrational_number_prime():
    """
    Finds the 6th prime in the concatenated digits of e and verifies its properties.
    """
    # Set a high precision to ensure we find the prime number, which is 1248 digits long.
    decimal.getcontext().prec = 1300

    # Get the digits of e, starting with '2'
    e_digits = str(decimal.Decimal(1).exp()).replace('.', '')

    primes_found_count = 0
    required_ending = 521023
    
    # Iterate through the digits of e, forming numbers to test for primality
    for i in range(1, len(e_digits)):
        current_num_str = e_digits[:i]
        
        # A small optimization: most large primes don't end in an even digit or 5
        if len(current_num_str) > 1 and current_num_str[-1] in '024568':
            continue

        current_num = int(current_num_str)
        
        if is_prime(current_num):
            primes_found_count += 1
            if primes_found_count == 6:
                sixth_prime = current_num
                break
    else:
        print("Could not find the 6th prime with the current precision.")
        return

    # Verify that the last 6 digits match the requirement
    if int(str(sixth_prime)[-6:]) == required_ending:
        print("The irrational number is Euler's number, e.")
        
        # As requested, construct and print an equation showing the full prime.
        # Equation form: (Prime_Stem * 1,000,000) + Last_6_Digits = Full_Prime
        prime_stem = sixth_prime // 1000000
        last_digits = required_ending
        
        print("\nThe equation demonstrating this is:")
        print(f"{prime_stem} * 1000000 + {last_digits} = {sixth_prime}")
    else:
        print(f"The 6th prime found ({sixth_prime}) does not end in {required_ending}.")


find_irrational_number_prime()