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

# Count prime numbers between 1825 and 2029
prime_count = 0
for number in range(1825, 2030):
    if is_prime(number):
        prime_count += 1

print(prime_count)