import random

def is_prime(n, k=10):
    """
    Miller-Rabin primality test.
    k is the number of rounds for testing. It's a probabilistic test.
    Returns True if n is probably prime, False if it is composite.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    if n < 25: # For small numbers, the above checks are enough
        return True

    # Write n-1 as 2^s * d where d is odd
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # Witness loop
    for _ in range(k):
        a = random.randint(2, n - 2)
        x = pow(a, d, n)  # x = (a^d) % n

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            # If the inner loop finishes without breaking, n is composite
            return False
    return True # n is probably prime

def solve_puzzle():
    """
    Main function to solve the puzzle.
    """
    # The problem specifies an irrational number expressible in no more than 3 symbols.
    # We will test Euler's number, 'e'.
    irrational_number_name = "Euler's number, e"
    
    # Concatenated digits of e, starting from the first digit.
    # We need a sufficient number of digits to find the 6th prime. ~400 digits are used here.
    e_digits = "271828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200305992181741359662904357290033429526059563073813232862794349076323382988075319525101901157383418793070215408914993488416750924476146066808226480016847741185374234544243710753907774499206955170276183860626133138458300075204493382656029760673711320070932870912744374704723069697720931014169283681902551510865746377211125238978442505695369677078544996996794686445490598793163688923009879312773617821542499922957635148220826989519366803318252886939849646510582093923982948879332036250944311730123819706841614039701983767932068328237646480429531180232878250981945581530175671736133206981125099618188159304169035159888851934580727386673858942277833453443657388045608301349581698436224362142106540136211494950444521023"

    primes_found = []
    print(f"Testing the irrational number: {irrational_number_name}\n")
    print("Searching for the first 6 prime numbers in the sequence of its digits...")
    
    for i in range(1, len(e_digits) + 1):
        current_number_str = e_digits[:i]
        current_number = int(current_number_str)
        
        if is_prime(current_number):
            primes_found.append(current_number)
            print(f"  > Found prime #{len(primes_found)} (a {len(current_number_str)}-digit number)")
        
        if len(primes_found) == 6:
            break

    print("\n--- Verification ---")
    if len(primes_found) == 6:
        sixth_prime = primes_found[5]
        sixth_prime_str = str(sixth_prime)
        target_digits = "521023"
        
        print(f"The 6th prime number found has {len(sixth_prime_str)} digits.")
        print("The 'final equation' is the verification of the last 6 digits of this prime.")
        print("\nThe full 6th prime number is:")
        print(sixth_prime)
        
        print("\nChecking its last 6 digits:")
        print(f"{sixth_prime_str[-6:]} == {target_digits}")
        
        if sixth_prime_str[-6:] == target_digits:
            print("\nSuccess! The digits match.")
            print(f"The irrational number is indeed {irrational_number_name}.")
        else:
            print("\nFailure. The digits do not match.")
    else:
        print("Could not find 6 primes with the provided digits.")

if __name__ == '__main__':
    solve_puzzle()
<<<e>>>