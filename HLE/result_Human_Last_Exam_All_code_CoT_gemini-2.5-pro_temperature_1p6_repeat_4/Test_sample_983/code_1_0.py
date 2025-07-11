import math

def get_sum_of_divisors(n):
    """Calculates the sum of all divisors of a number n (sigma(n))."""
    if n <= 0:
        return 0
    sum_div = 0
    # Iterating up to the square root is an efficient way to find divisors
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            # i is a divisor
            sum_div += i
            # n/i is also a divisor
            if i * i != n:
                sum_div += n // i
    return sum_div

def is_zumkeller(n):
    """
    Checks if a number n is a Zumkeller number based on its properties.
    """
    sigma_n = get_sum_of_divisors(n)
    
    # Condition 1: The number must be perfect or abundant.
    is_abundant_or_perfect = (sigma_n >= 2 * n)
    
    # Condition 2: The number and the sum of its divisors must have the same parity.
    has_same_parity = (sigma_n % 2 == n % 2)
    
    return is_abundant_or_perfect and has_same_parity

def find_next_in_sequence():
    """
    Finds the next Zumkeller number after the last term in the given sequence.
    """
    last_term_in_sequence = 2352
    
    print(f"The sequence is composed of Zumkeller numbers. We are looking for the next one after {last_term_in_sequence}.")
    
    candidate_number = last_term_in_sequence + 1
    
    while True:
        if is_zumkeller(candidate_number):
            next_number = candidate_number
            break
        candidate_number += 1
        
    # As requested, printing the "equation" part by showing the verification steps.
    sigma_next = get_sum_of_divisors(next_number)
    print(f"\nFound the next number: {next_number}")
    print("\n--- Verification ---")
    print(f"For n = {next_number}:")
    print(f"Sum of divisors, sigma({next_number}) = {sigma_next}")
    print(f"Twice the number, 2 * {next_number} = {2 * next_number}")
    
    print("\nChecking Zumkeller conditions:")
    print(f"1. Is sigma(n) >= 2*n?  =>  Is {sigma_next} >= {2 * next_number}?  Result: {sigma_next >= 2 * next_number}")
    print(f"2. Is parity of n and sigma(n) the same? => ({next_number} % 2) == ({sigma_next} % 2)? Result: {next_number % 2 == sigma_next % 2}")
    
    print(f"\nSince both conditions are met, the next number in the sequence is {next_number}.")

if __name__ == '__main__':
    find_next_in_sequence()