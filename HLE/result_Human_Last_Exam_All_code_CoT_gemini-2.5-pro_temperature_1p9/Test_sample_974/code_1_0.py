import math

def get_prime_divisors(n):
    """
    Returns a list of unique prime divisors of a number n.
    """
    factors = set()
    # Check for factor 2
    if n % 2 == 0:
        factors.add(2)
        while n % 2 == 0:
            n = n // 2
    # Check for odd factors
    i = 3
    limit = int(math.sqrt(n))
    while i <= limit:
        if n % i == 0:
            factors.add(i)
            while n % i == 0:
                n = n // i
            limit = int(math.sqrt(n)) # Update limit
        else:
            i += 2
    if n > 1:
        factors.add(n)
    return sorted(list(factors))

def solve_and_print():
    """
    Solves the problem as described and prints the detailed steps and final answer.
    """
    q = 12740347

    print("Step-by-step solution:")
    print(f"Let q = {q}.")
    print("The problem asks for all 'primes divisors p' for which the number of elements of order p is equal in PSL(3, q^2) and PSL(4, q).")
    print("This is interpreted as finding prime divisors p of the number q that satisfy the condition.")
    
    print(f"\nStep 1: Find all prime divisors of q = {q}.")
    prime_divs_of_q = get_prime_divisors(q)
    
    if len(prime_divs_of_q) == 1 and prime_divs_of_q[0] == q:
        print(f"The only prime divisor of {q} is {q} itself, which means q is a prime number.")
    else:
        print(f"The prime divisors of {q} are: {prime_divs_of_q}")

    print("\nStep 2: For each prime divisor p, check the condition.")
    solutions = []
    for p in prime_divs_of_q:
        print(f"\n--- Checking for p = {p} ---")
        # Since q is prime, the only prime divisor p is q itself.
        # The characteristic of the fields GF(q^2) and GF(q) is q.
        # We are counting elements of order p=q, which are unipotent elements.

        print("The number of elements of order p=q is the number of non-identity unipotent elements.")
        print("The number of unipotent elements in SL(n, k) is k^(n*(n-1)).")
        
        # For PSL(3, q^2)
        n1 = 3
        power_1 = n1 * (n1 - 1) * 2
        
        # For PSL(4, q)
        n2 = 4
        power_2 = n2 * (n2 - 1)

        print("\nFor PSL(3, q^2):")
        print(f"Here, n = {n1} and the field is GF(q^2), so the base k is q^2.")
        print(f"Number of elements of order p = (q^2)^({n1}*({n1-1})) - 1 = q^({2*n1*(n1-1)}) - 1 = q^{power_1} - 1.")

        print("\nFor PSL(4, q):")
        print(f"Here, n = {n2} and the field is GF(q), so the base k is q.")
        print(f"Number of elements of order p = q^({n2}*({n2-1})) - 1 = q^{power_2} - 1.")

        print("\nComparison:")
        print(f"The number of elements for PSL(3, q^2) is q^{power_1} - 1.")
        print(f"The number of elements for PSL(4, q) is q^{power_2} - 1.")

        if power_1 == power_2:
            print(f"The exponents are equal ({power_1} == {power_2}), so the counts are equal.")
            print(f"Therefore, p = {p} is a solution.")
            solutions.append(p)
        else:
            print(f"The exponents are not equal ({power_1} != {power_2}), so p = {p} is not a solution.")

    print("\n" + "="*20 + " Final Answer " + "="*20)
    if not solutions:
        print("No prime divisors of q satisfy the condition.")
    else:
        print("The prime divisor p which satisfies the condition is:")
        result = solutions[0]
        print(result)
            
        print("\nThe final equality of the number of elements is expressed by the equation:")
        n1, n2 = 3, 4
        # Left side for PSL(3, q^2)
        print(f"({result}^2)^({n1} * ({n1} - 1)) - 1 = ", end="")
        # Right side for PSL(4, q)
        print(f"{result}^({n2} * ({n2} - 1)) - 1")

        # Simplified form
        print(f"Which simplifies to:")
        print(f"{result}^{power_1} - 1 = {result}^{power_2} - 1")

solve_and_print()
<<<12740347>>>