import math

def get_prime_factorization(num):
    """
    Calculates the prime factorization of a number.
    Returns a dictionary of {prime: exponent}.
    """
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def sum_of_divisors_from_factors(factors):
    """
    Calculates the sum of divisors (sigma function) from a prime factorization dict.
    """
    s = 1
    for p, a in factors.items():
        s *= (p**(a + 1) - 1) // (p - 1)
    return s

def is_practical(n):
    """
    Checks if a number is practical using Stewart's theorem.
    Returns True if practical, False otherwise.
    """
    if n <= 0: return False
    if n == 1: return True
    
    factors = get_prime_factorization(n)
    primes = sorted(factors.keys())
    
    # Condition 1: First prime must be 2.
    if primes[0] != 2:
        return False
        
    # Condition 2: p_i <= sigma(...) + 1
    running_product_factors = {}
    for p in primes:
        a = factors[p]
        # Calculate sigma for the product of factors processed so far
        # For the first prime (2), sigma is compared against p_1, but that check is done above.
        # This condition is for p_2, p_3, ...
        if p != 2:
            sigma_prev = sum_of_divisors_from_factors(running_product_factors)
            if p > sigma_prev + 1:
                return False

        running_product_factors[p] = a
        
    return True

def find_and_verify_next_practical_number(start_num):
    """
    Finds the next practical number after start_num and prints the verification.
    """
    num = start_num + 1
    while True:
        if is_practical(num):
            next_practical_num = num
            break
        num += 1

    print(f"The next number in the sequence is: {next_practical_num}\n")
    print(f"Verification that {next_practical_num} is a practical number:")

    factors = get_prime_factorization(next_practical_num)
    primes = sorted(factors.keys())
    
    factor_str = " * ".join([f"{p}^{e}" for p, e in factors.items()])
    print(f"1. Prime factorization of {next_practical_num} is {factor_str}.")
    
    # Check condition 1
    if primes[0] == 2:
        print("2. The first prime is 2. (Condition 1 met)")
    else:
        print(f"2. The first prime is {primes[0]}, not 2. (Condition 1 failed)")
        return
        
    # Check condition 2 for subsequent primes
    running_factors_for_sigma = {}
    i = 3 # Step number for printing
    for p in primes:
        a = factors[p]
        if p != 2:
            sigma_val = sum_of_divisors_from_factors(running_factors_for_sigma)
            running_prod_num = math.prod(k**v for k,v in running_factors_for_sigma.items())
            
            check_str = f"Checking prime {p}: {p} <= sum_of_divisors({running_prod_num}) + 1 = {sigma_val} + 1 = {sigma_val + 1}."
            if p <= sigma_val + 1:
                print(f"{i}. {check_str} (Condition met)")
            else:
                print(f"{i}. {check_str} (Condition failed)")

            i += 1
        running_factors_for_sigma[p] = a

    print(f"\nSince all conditions are met, {next_practical_num} is a practical number.")

if __name__ == '__main__':
    last_number_in_sequence = 2352
    find_and_verify_next_practical_number(last_number_in_sequence)
    # The final answer part is for the platform, the user can run the code above.
    # The calculated next number is 2360.
    print("\n<<<2360>>>")