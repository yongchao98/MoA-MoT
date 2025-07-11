import math

def find_prime_divisors(n):
    """
    Finds all unique prime divisors of a given integer n.
    Returns a sorted list of the prime factors.
    """
    factors = set()
    # Check for divisibility by 2
    if n % 2 == 0:
        factors.add(2)
        while n % 2 == 0:
            n //= 2
    # Check for odd divisors from 3 upwards
    d = 3
    # We only need to check up to the square root of the number
    limit = int(math.sqrt(n))
    while d <= limit:
        if n % d == 0:
            factors.add(d)
            while n % d == 0:
                n //= d
            # Update the limit since n has been reduced
            limit = int(math.sqrt(n))
        d += 2
    # If n is still greater than 1 at this point, it must be a prime factor
    if n > 1:
        factors.add(n)
    return sorted(list(factors))

def solve_and_print():
    """
    Calculates the solution and prints the steps as requested.
    """
    q = 12740347

    print(f"The given number is q = {q}.")
    print("The goal is to find all prime divisors p of q such that the number of elements of order p in PSL(3, q^2) and PSL(4, q) are equal.")
    print("\nThe number of elements of order p (where p is the field characteristic) in PSL(n, k) is k^(n*(n-1)) - 1.")

    print("\nFor the group PSL(3, q^2):")
    n1 = 3
    print(f"The parameter n is {n1}.")
    print(f"The field size k is q^2.")
    print(f"Number of elements = (q^2)^(n*(n-1)) - 1 = (q^2)^({n1}*({n1}-1)) - 1 = (q^2)^({n1*2}) - 1 = q^12 - 1.")

    print("\nFor the group PSL(4, q):")
    n2 = 4
    print(f"The parameter n is {n2}.")
    print(f"The field size k is q.")
    print(f"Number of elements = q^(n*(n-1)) - 1 = q^({n2}*({n2}-1)) - 1 = q^({n2*3}) - 1 = q^12 - 1.")

    print("\nSince both formulas evaluate to q^12 - 1, the condition holds for all prime divisors of q.")
    print("We now find all prime divisors of q.")

    prime_divisors = find_prime_divisors(q)

    print("\nThe list of prime divisors p is:")
    for p in prime_divisors:
        print(p)

solve_and_print()
<<<12740347>>>