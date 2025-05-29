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

def count_primes_in_range(start, end):
    prime_count = 0
    for num in range(start, end + 1):
        if is_prime(num):
            prime_count += 1
    return prime_count

# Define the range
start = 8667
end = 9398

# Count the prime numbers in the range
prime_count = count_primes_in_range(start, end)
print(prime_count)