def is_prime(n):
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

def largest_prime_or_number(numbers):
    primes = [num for num in numbers if is_prime(num)]
    if primes:
        return max(primes)
    else:
        return max(numbers)

input_list = [49, 48, 93]
output = largest_prime_or_number(input_list)
print(output)