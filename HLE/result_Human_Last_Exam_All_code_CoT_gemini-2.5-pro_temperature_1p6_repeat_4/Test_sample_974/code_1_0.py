import math

def get_prime_divisors(n):
    """
    Returns a sorted list of unique prime factors of a natural number n.
    """
    factors = set()
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.add(d)
            while temp_n % d == 0:
                temp_n //= d
        d += 1
    if temp_n > 1:
       factors.add(temp_n)
    return sorted(list(factors))

def solve():
    """
    Finds all prime divisors p of q=12740347 for which the number of elements
    of order p in PSL(3, q^2) and PSL(4, q) are equal.
    """
    q = 12740347

    # Step 1: Find all prime divisors of q. For F_q to be a field,
    # q must be a prime power.
    print(f"Let q = {q}.")
    print("Step 1: Find the prime divisors of q.")
    prime_divisors = get_prime_divisors(q)

    # We find that q is a prime number.
    if len(prime_divisors) == 1 and prime_divisors[0] == q:
        print(f"The number q = {q} is prime.")
        print(f"The only prime divisor of q is p = {q}.")
        # The characteristic of the field F_q is q.
    else:
        # This case is not applicable here but included for completeness.
        print(f"The prime divisors of q are: {prime_divisors}")
    
    # We only need to check for p = q.
    p = q
    print(f"\nThus, we need to check the condition for the prime p = {p}.")
    print("This prime p is the characteristic of the finite fields F_q and F_{q^2}.")

    # Step 2: Use the formula for the number of elements of order p (the characteristic)
    # in PSL(n, F_r), which is N_p = r^(n*(n-1)/2) - 1.
    print("\nStep 2: Calculate the number of elements of order p for each group.")

    # For PSL(3, q^2): n=3, r = q^2
    n1 = 3
    print(f"For PSL(3, q^2), n=3 and r = q^2 = {q}^2.")
    print(f"The number of elements of order {p} is r^(n(n-1)/2) - 1 = (q^2)^({n1}({n1}-1)/2) - 1")
    print(f"= (q^2)^3 - 1 = q^6 - 1.")

    # For PSL(4, q): n=4, r = q
    n2 = 4
    print(f"\nFor PSL(4, q), n=4 and r = q = {q}.")
    print(f"The number of elements of order {p} is r^(n(n-1)/2) - 1 = q^({n2}({n2}-1)/2) - 1")
    print(f"= q^6 - 1.")

    # Step 3: Compare the two quantities.
    print("\nStep 3: Compare the two results.")
    print(f"The condition is that the number of elements are equal:")
    print(f"Number in PSL(3, {q}^2) = Number in PSL(4, {q})")
    print(f"{q}^6 - 1 = {q}^6 - 1")
    print("This equality is true for any q > 1.")

    # Step 4: Final Conclusion
    print("\nStep 4: Conclusion.")
    print(f"The condition holds for the prime divisor p = {q}.")
    print("Therefore, the list of all prime divisors of q satisfying the condition is:")
    print(prime_divisors)

if __name__ == '__main__':
    solve()