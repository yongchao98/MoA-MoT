import math

def get_prime_divisors(n):
    """
    Returns a list of unique prime divisors of a number n.
    """
    factors = set()
    d = 2
    temp_n = n
    # Test divisibility up to sqrt(n)
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.add(d)
            while temp_n % d == 0:
                temp_n //= d
        d += 1
    # The remaining number is also a prime factor
    if temp_n > 1:
        factors.add(temp_n)
    return list(factors)

def solve_and_print():
    """
    Solves the problem and prints the derivation and final answer.
    """
    q = 12740347

    print(f"Given the number q = {q}.")
    
    # Step 1: Find the prime divisors of q.
    prime_divisors = get_prime_divisors(q)
    print(f"The prime divisors of q are: {prime_divisors}")
    print("-" * 30)

    # Since q is prime, its only prime divisor is q itself.
    # We only need to test p = q.
    p = prime_divisors[0]
    
    print(f"We need to check if the number of elements of order p = {p} are equal in both groups.")
    
    # Step 2: Explain the formula for the number of elements.
    print("\nThe elements of order p (where p is the characteristic of the field) are unipotent.")
    print("For a group PSL(n, k) with k being a field of characteristic p:")
    print("If p > n, all non-identity unipotent elements have order exactly p.")
    print("The number of such elements is k^(n*(n-1)/2) - 1.")
    print("-" * 30)

    # Step 3: Calculate the number for the first group, PSL(3, q^2).
    n1 = 3
    print(f"For PSL({n1}, q^2):")
    print(f"  - The field is GF(q^2) = GF({q}^2), its characteristic is p = {p}.")
    print(f"  - The dimension is n = {n1}.")
    print(f"  - The condition p > n ({p} > {n1}) is met.")
    print(f"  - The number of elements of order {p} is:")
    # Using format to show the equation as requested
    print(f"    (q^2)^(n*(n-1)/2) - 1 = ({p}^2)^({n1}*({n1}-1)/2) - 1 = ({p}^2)^3 - 1 = {p}^6 - 1")

    # Step 4: Calculate the number for the second group, PSL(4, q).
    n2 = 4
    print(f"\nFor PSL({n2}, q):")
    print(f"  - The field is GF(q) = GF({p}), its characteristic is p = {p}.")
    print(f"  - The dimension is n = {n2}.")
    print(f"  - The condition p > n ({p} > {n2}) is met.")
    print(f"  - The number of elements of order {p} is:")
    # Using format to show the equation as requested
    print(f"    q^(n*(n-1)/2) - 1 = {p}^({n2}*({n2}-1)/2) - 1 = {p}^6 - 1")

    # Step 5: Compare the results and conclude.
    print("-" * 30)
    print("Comparing the two results, we see that both are equal to q^6 - 1.")
    print(f"Therefore, the prime p = {p} satisfies the condition.")

    print("\nThe list of all such prime divisors is:")
    print(p)

solve_and_print()