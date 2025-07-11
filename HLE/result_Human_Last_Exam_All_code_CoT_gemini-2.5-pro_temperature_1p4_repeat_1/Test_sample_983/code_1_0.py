import collections

def get_prime_factorization(num):
    """
    Returns a dictionary of prime factors of num and their counts.
    For example, get_prime_factorization(12) returns {2: 2, 3: 1}.
    """
    factors = collections.defaultdict(int)
    d = 2
    temp = num
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def is_sphenic(n):
    """
    Checks if a number is sphenic.
    A sphenic number is a product of three distinct prime numbers.
    """
    if n <= 1:
        return False
    factors = get_prime_factorization(n)
    # A sphenic number must be a product of exactly 3 primes.
    if len(factors) != 3:
        return False
    # Each prime factor must appear exactly once (distinct primes).
    if any(count != 1 for count in factors.values()):
        return False
    return True

def sum_digits(n):
    """Calculates the sum of the digits of a number."""
    s = 0
    for digit in str(n):
        s += int(digit)
    return s

def find_next_number_in_sequence(last_number):
    """
    Finds the next 'mod-sphenic' number after last_number.
    A number n is 'mod-sphenic' if n + sum_digits(n) is a sphenic number.
    """
    current_num = last_number + 1
    while True:
        num_to_check = current_num + sum_digits(current_num)
        if is_sphenic(num_to_check):
            return current_num
        current_num += 1

# The last number in the sequence is 2352.
last_number_in_sequence = 2352
next_number = find_next_number_in_sequence(last_number_in_sequence)

# Now, we format the output as requested.
s_digits = sum_digits(next_number)
resulting_sphenic_number = next_number + s_digits
factors = get_prime_factorization(resulting_sphenic_number)
p1, p2, p3 = sorted(factors.keys())

# Print the explanation of the result.
print(f"The rule is that for a number 'n' in the sequence, 'n + sum of its digits' must be a sphenic number (a product of three distinct primes).")
print(f"The last number provided is {last_number_in_sequence}.")
print(f"Searching for the next number...\n")
print(f"The next number in the sequence is {next_number}.")
print(f"This is because:")
# Display the sum of digits calculation
digit_sum_str = "+".join(str(next_number))
print(f"Sum of digits of {next_number} = {'+'.join(list(str(next_number)))} = {s_digits}")
# Display the equation that results in a sphenic number
print(f"{next_number} + {s_digits} = {resulting_sphenic_number}")
# Display the prime factorization that proves the number is sphenic
print(f"And {resulting_sphenic_number} is a sphenic number, as it is the product of three distinct primes: {p1} * {p2} * {p3}.")

# The final answer in the required format
print(f"\n<<<2359>>>")