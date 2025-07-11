import math

def get_prime_factorization(n):
    """Returns a dictionary of prime factors and their counts."""
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def phi(n):
    """Computes Euler's totient function phi(n)."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    result = n
    for p in factors:
        result -= result // p
    return result

def get_divisors(n):
    """Returns a sorted list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve():
    """
    Solves the problem of finding the number of primitive Dirichlet characters
    of a given conductor and order.
    """
    d = 53599
    k = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order k = {k}.")
    print("=" * 70)

    # Step 1: Factorize d
    prime_factors = get_prime_factorization(d)
    primes = sorted(list(prime_factors.keys()))
    print(f"Step 1: Analyze the conductor d = {d}")
    print(f"The prime factorization of d is: {d} = {' * '.join(map(str, primes))}")
    print("Since d is square-free, a character chi mod d is primitive if and only if its")
    print("component chi_p for each prime factor p is itself primitive (i.e., non-trivial, order > 1).\n")

    # Step 2: Explain the counting method using Inclusion-Exclusion
    print("Step 2: Method for counting characters of a specific order")
    print("The order of a character chi = product(chi_p) is lcm(order(chi_p)). We need this lcm to be 6.")
    print("Let A(m) be the number of primitive character combinations where the order of each component divides m.")
    print("By the principle of inclusion-exclusion, the number of characters with order exactly 6 is:")
    print("Result = A(6) - A(3) - A(2) + A(1)\n")
    
    print("Step 3: Calculate the terms A(m)")
    # Check if for all p, (p-1) is divisible by all j that divide k.
    # For this problem, k=6. The primes are 7, 13, 19, 31.
    # p-1 values are 6, 12, 18, 30. All are divisible by 6.
    # So any divisor of 6 also divides any p-1.

    terms = {}
    divisors_of_k = get_divisors(k)

    for m in divisors_of_k:
        # Get divisors of m, must be > 1 for primitive characters.
        orders_to_consider = [j for j in get_divisors(m) if j > 1]
        
        # Calculate number of character choices for a single prime factor.
        # This is sum of phi(j) for all relevant orders j.
        num_choices_per_prime = sum(phi(j) for j in orders_to_consider)
        
        # Calculate total combinations A(m) for all prime factors.
        num_primes = len(primes)
        a_m = num_choices_per_prime ** num_primes
        terms[m] = a_m

        print(f"Calculating A({m}):")
        if not orders_to_consider:
            print("  - Component character orders must be > 1 and divide {m}. No such orders exist.")
        else:
            phi_calcs = ' + '.join([f'phi({j})' for j in orders_to_consider])
            phi_vals = ' + '.join([str(phi(j)) for j in orders_to_consider])
            print(f"  - Component character orders must be in {orders_to_consider}.")
            print(f"  - Number of choices per prime = {phi_calcs} = {phi_vals} = {num_choices_per_prime}.")
        
        print(f"  - With {num_primes} prime factors, A({m}) = {num_choices_per_prime}^{num_primes} = {a_m}.\n")
        
    print("=" * 70)
    print("Step 4: Final Calculation")
    
    A6 = terms.get(6, 0)
    A3 = terms.get(3, 0)
    A2 = terms.get(2, 0)
    A1 = terms.get(1, 0)

    final_result = A6 - A3 - A2 + A1
    
    print(f"The number of primitive characters of order 6 is A(6) - A(3) - A(2) + A(1).")
    print(f"Result = {A6} - {A3} - {A2} + {A1}")
    print(f"Result = {A6 - A3 - A2 + A1}")
    print(f"\nThus, there are {final_result} primitive Dirichlet characters of conductor {d} and order {k}.")
    
    return final_result

# Run the solver and store the final answer
final_answer = solve()

# The final answer is wrapped according to the required format
# print(f"\n<<<608>>>")