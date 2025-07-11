import decimal
import random

def is_prime(n, k=5):
    """
    Miller-Rabin primality test. A probabilistic test that is very accurate for large numbers.
    k is the number of rounds of testing to perform.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False

    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1

    for _ in range(k):
        a = random.randrange(2, n - 2)
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

def solve_riddle():
    """
    Finds the 6th prime in the sequence of concatenated digits of 'e'
    and verifies it ends with 521023.
    """
    # The 6th prime is 986 digits long, so we need at least that much precision.
    # 1000 should be safe.
    precision = 1000
    decimal.getcontext().prec = precision + 2  # Add a buffer for precision

    print("The irrational number is 'e'. Its digits start 2.71828...")
    print(f"Calculating the first {precision} digits of 'e' and searching for primes...")
    print("-" * 30)

    # Get the digits of 'e' as a string, removing the decimal point.
    e_digits = str(decimal.Decimal(1).exp()).replace('.', '')

    primes_found = []
    
    # Iterate through prefixes of the digit string
    for i in range(1, len(e_digits) + 1):
        num_str = e_digits[:i]
        
        # Small optimization: skip numbers ending in 0, 2, 4, 5, 6, 8
        if len(num_str) > 1 and int(num_str[-1]) in [0, 2, 4, 5, 6, 8]:
            continue

        current_num = int(num_str)

        if is_prime(current_num):
            primes_found.append(current_num)
            print(f"Found prime #{len(primes_found)} (with {len(num_str)} digits): {current_num}")
            if len(primes_found) == 6:
                break
    
    print("-" * 30)

    if len(primes_found) == 6:
        sixth_prime = primes_found[5]
        last_six_digits = str(sixth_prime)[-6:]
        
        print("Final Answer Verification:")
        print(f"The 6th prime number in the sequence is:")
        print(sixth_prime)
        
        print("\nThe final equation is checking the last 6 digits:")
        print(f"{sixth_prime} % 1000000 = {last_six_digits}")

        if last_six_digits == "521023":
            print("\nThe last 6 digits match the number 521023 from the riddle.")
        else:
            print("\nVerification failed. The last 6 digits do not match.")
    else:
        print("Could not find 6 primes. Please increase the precision.")

if __name__ == "__main__":
    solve_riddle()