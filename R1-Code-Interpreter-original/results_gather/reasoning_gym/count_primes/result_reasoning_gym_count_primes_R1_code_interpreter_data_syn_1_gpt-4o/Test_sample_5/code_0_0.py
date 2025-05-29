def is_prime(n):
    """Check if a number is prime."""
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

# List of numbers to check
numbers = range(9958, 9964)

# Count prime numbers in the list
prime_count = sum(1 for number in numbers if is_prime(number))

print(prime_count)