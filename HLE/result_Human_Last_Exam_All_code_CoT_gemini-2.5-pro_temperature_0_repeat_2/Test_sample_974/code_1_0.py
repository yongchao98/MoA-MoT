import math

def find_prime_divisors(n):
    """
    Finds all unique prime divisors of a given integer n.
    """
    divisors = set()
    # Check for divisibility by 2
    if n % 2 == 0:
        divisors.add(2)
        while n % 2 == 0:
            n //= 2
    # Check for odd divisors from 3 up to sqrt(n)
    d = 3
    limit = int(math.sqrt(n))
    while d <= limit:
        if n % d == 0:
            divisors.add(d)
            while n % d == 0:
                n //= d
            # Update the limit as n changes
            limit = int(math.sqrt(n))
        d += 2
    # If n is still greater than 1, it must be a prime factor itself
    if n > 1:
        divisors.add(n)
    return sorted(list(divisors))

# Given value of q
q = 12740347

# For PSL(3, q^2): n=3, k=q^2
n1 = 3
exp1 = n1 * (n1 - 1)
final_exp1 = 2 * exp1

# For PSL(4, q): n=4, k=q
n2 = 4
exp2 = n2 * (n2 - 1)

print(f"Let q = {q}.")
print("The number of elements of order p (where p is a prime divisor of q) is determined by the group structure.")
print("\nFor PSL(3, q^2):")
print(f"The formula is k^(n(n-1)) - 1, with n = {n1} and k = q^2.")
print(f"Number of elements = (q^2)^({n1}({n1}-1)) - 1 = (q^2)^{exp1} - 1 = q^{final_exp1} - 1.")

print("\nFor PSL(4, q):")
print(f"The formula is k^(n(n-1)) - 1, with n = {n2} and k = q.")
print(f"Number of elements = q^({n2}({n2}-1)) - 1 = q^{exp2} - 1.")

print(f"\nSince q^{final_exp1} - 1 = q^{exp2} - 1, the number of elements is equal.")
print("This holds for any prime divisor p of q.")

# Find and print the prime divisors of q
prime_divisors_of_q = find_prime_divisors(q)

print(f"\nThe prime divisors of q = {q} are:")
for p in prime_divisors_of_q:
    print(p)