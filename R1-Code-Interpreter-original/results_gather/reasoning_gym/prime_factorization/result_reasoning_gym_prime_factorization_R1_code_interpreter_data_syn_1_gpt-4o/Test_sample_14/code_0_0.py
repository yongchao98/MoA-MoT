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

def prime_factors(n):
    factors = []
    # Check for number of 2s
    while n % 2 == 0:
        factors.append(2)
        n = n // 2
    # Check for odd factors
    for i in range(3, int(n**0.5) + 1, 2):
        while n % i == 0:
            factors.append(i)
            n = n // i
    # If n is a prime number greater than 2
    if n > 2:
        factors.append(n)
    return factors

number = 911
if is_prime(number):
    print(f"{number} is a prime number.")
else:
    factors = prime_factors(number)
    print(" × ".join(map(str, factors)))
