import math

def calculate_n_for_construction(prime_exponent_ranges):
    """
    Calculates the maximum integer N for a given construction.
    prime_exponent_ranges is a list of tuples (prime, [exponents]).
    """
    max_val = 0
    
    def generate_numbers(primes, ranges, current_num=1):
        nonlocal max_val
        if not primes:
            if current_num > max_val:
                max_val = current_num
            return

        prime = primes[0]
        remaining_primes = primes[1:]
        exponent_range = ranges[0]
        remaining_ranges = ranges[1:]
        
        for exp in exponent_range:
            generate_numbers(remaining_primes, remaining_ranges, current_num * (prime**exp))

    primes = [item[0] for item in prime_exponent_ranges]
    ranges = [item[1] for item in prime_exponent_ranges]
    
    generate_numbers(primes, ranges)
    return max_val

# Case 1: 16 -> One prime {2}, exponents {0,...,15}
n1 = calculate_n_for_construction([(2, range(16))])
print(f"Construction 1 (1 prime, 16 exponents): 2^15 = {n1}")

# Case 2: 4x4 -> Two primes {2,3}, exponents {0,1,2,3} for each
n2 = calculate_n_for_construction([(2, range(4)), (3, range(4))])
print(f"Construction 2 (2 primes, 4x4 exponents): max(2^a * 3^b) = 2^3 * 3^3 = {n2}")

# Case 3: 2x2x2x2 -> Four primes {2,3,5,7}, exponents {0,1} for each
n3 = calculate_n_for_construction([(2, range(2)), (3, range(2)), (5, range(2)), (7, range(2))])
print(f"Construction 3 (4 primes, 2x2x2x2 exponents): 2*3*5*7 = {n3}")

# Case 4: 4x2x2 -> Three primes {2,3,5}, exponents {0,1,2,3}, {0,1}, {0,1}
n4 = calculate_n_for_construction([(2, range(4)), (3, range(2)), (5, range(2))])
print(f"Construction 4 (3 primes, 4x2x2 exponents): max(2^a * 3^b * 5^c) = 2^3 * 3^1 * 5^1 = {n4}")

min_n = min(n1, n2, n3, n4)
print(f"\nThe smallest N among these constructions is: {min_n}")