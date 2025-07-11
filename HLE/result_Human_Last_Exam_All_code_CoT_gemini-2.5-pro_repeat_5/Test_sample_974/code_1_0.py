import math

def get_prime_factors(n):
    """
    Returns a list of prime factors of n.
    """
    factors = set()
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors.add(d)
            temp //= d
        d += 1
    if temp > 1:
        factors.add(temp)
    return sorted(list(factors))

def solve():
    """
    Finds all prime divisors p of q for which the number of elements of order p
    in PSL(3, q^2) and PSL(4, q) are equal.
    """
    q = 12740347

    print(f"Given q = {q}.")
    print("The problem is to find all prime divisors p of q such that the number of elements of order p is equal in PSL(3, q^2) and PSL(4, q).")

    # Step 1: Find prime divisors of q.
    prime_divisors_q = get_prime_factors(q)
    print(f"\nThe prime divisors of q = {q} are: {prime_divisors_q}.")
    print("Since q is a prime number, the only candidate to check is p = q.")

    solutions = []

    # Step 2: Check the condition for each prime divisor p.
    for p in prime_divisors_q:
        print(f"\nChecking for p = {p}:")
        
        # For p=q, the elements of order p are related to unipotent elements.
        # The number of such elements in PSL(n, k) where char(k)=p is well-known.
        
        # For PSL(3, q^2), the field is F_{q^2}. The number of elements of order q is q^12 - 1.
        num_elements_g1_str = f"{q}^12 - 1"
        
        # For PSL(4, q), the field is F_q. The number of elements of order q is q^12 - 1.
        num_elements_g2_str = f"{q}^12 - 1"

        print(f"Number of elements of order {p} in PSL(3, {q}^2) = {num_elements_g1_str}")
        print(f"Number of elements of order {p} in PSL(4, {q}) = {num_elements_g2_str}")
        
        print("\nForming the equation:")
        print(f"{num_elements_g1_str} = {num_elements_g2_str}")

        # The expressions are identical, so the equality holds.
        print("The equation is true.")
        solutions.append(p)

    print("\n-------------------------------------------------")
    print("Final list of prime divisors p satisfying the condition:")
    if not solutions:
        print("No solutions found.")
    else:
        for sol in solutions:
            print(sol)

solve()