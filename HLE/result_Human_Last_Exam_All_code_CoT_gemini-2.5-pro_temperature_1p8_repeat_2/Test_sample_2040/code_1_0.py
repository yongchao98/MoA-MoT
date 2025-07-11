import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        if n % i == 0:
            return False
    return True

def solve_task():
    """
    Finds the first 5 prime numbers, calculates their sum,
    and prints the full equation.
    """
    primes_to_find = 5
    primes = []
    num = 2
    while len(primes) < primes_to_find:
        if is_prime(num):
            primes.append(num)
        num += 1

    total_sum = sum(primes)

    # Format the equation string
    # The instructions require printing each number in the final equation.
    equation_str = " + ".join(map(str, primes))
    
    print(f"{equation_str} = {total_sum}")

solve_task()