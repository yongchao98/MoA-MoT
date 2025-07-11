import math

def prime_factorization(n):
    """Returns a dictionary of prime factors of n."""
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
    result = n
    for p in prime_factorization(n):
        result -= result // p
    return result

def solve():
    """
    Calculates the number of primitive Dirichlet characters of a given
    conductor and order.
    """
    d = 53599
    order_k = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order k = {d_factors}.")
    
    d_factors = prime_factorization(d)
    primes = list(d_factors.keys())
    print(f"Step 1: The prime factorization of d = {d} is { ' * '.join(map(str, primes)) }.")
    print("Since d is square-free, a character is primitive if and only if its components are non-principal (order > 1).")
    print(f"We need to find combinations of characters where the LCM of their orders is {order_k}.\n")

    print("Step 2: Use the inclusion-exclusion principle.")
    print("Result = (count where lcm divides 6) - (count where lcm divides 3) - (count where lcm divides 2)\n")

    # --- Calculate total where LCM divides 6 ---
    phi_2 = phi(2)
    phi_3 = phi(3)
    phi_6 = phi(6)
    num_per_prime_div6 = phi_2 + phi_3 + phi_6
    total_div6 = num_per_prime_div6 ** len(primes)
    
    print("Calculation for tuples where LCM divides 6:")
    print(f"  For each prime factor p, the number of primitive characters with order dividing 6 is:")
    print(f"  phi(2) + phi(3) + phi(6) = {phi_2} + {phi_3} + {phi_6} = {num_per_prime_div6}")
    print(f"  Total number of combinations = {num_per_prime_div6}^{len(primes)} = {total_div6}\n")

    # --- Calculate total where LCM divides 3 ---
    num_per_prime_div3 = phi(3)
    total_div3 = num_per_prime_div3 ** len(primes)
    
    print("Calculation for tuples where LCM divides 3:")
    print(f"  For each prime factor p, the number of primitive characters with order dividing 3 is:")
    print(f"  phi(3) = {phi_3}")
    print(f"  Total number of combinations = {num_per_prime_div3}^{len(primes)} = {total_div3}\n")

    # --- Calculate total where LCM divides 2 ---
    num_per_prime_div2 = phi(2)
    total_div2 = num_per_prime_div2 ** len(primes)
    
    print("Calculation for tuples where LCM divides 2:")
    print(f"  For each prime factor p, the number of primitive characters with order dividing 2 is:")
    print(f"  phi(2) = {phi_2}")
    print(f"  Total number of combinations = {num_per_prime_div2}^{len(primes)} = {total_div2}\n")

    # --- Final Calculation ---
    final_result = total_div6 - total_div3 - total_div2
    
    print("Step 3: Final Calculation")
    print("The total number of primitive characters of order 6 is:")
    print(f"   ({num_per_prime_div6}^{len(primes)}) - ({num_per_prime_div3}^{len(primes)}) - ({num_per_prime_div2}^{len(primes)})")
    print(f"=  {total_div6} - {total_div3} - {total_div2}")
    print(f"=  {final_result}")
    
    return final_result

if __name__ == '__main__':
    solve()
