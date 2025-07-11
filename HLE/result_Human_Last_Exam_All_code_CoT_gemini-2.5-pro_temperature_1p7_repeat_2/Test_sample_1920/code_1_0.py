import math

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary {prime: exponent}.
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def get_divisors(n):
    """
    Returns a sorted list of all divisors of n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def mobius(n):
    """
    Calculates the Mobius function mu(n).
    """
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    if any(exp > 1 for exp in factors.values()):
        return 0
    
    return -1 if len(factors) % 2 != 0 else 1

def main():
    """
    Calculates the number of primitive Dirichlet characters of a given
    conductor d and order k.
    """
    d = 53599
    k = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order k = {k}.")
    
    # Step 1: Factorize the conductor d.
    d_factors = get_prime_factorization(d)
    prime_factors = sorted(d_factors.keys())
    print(f"\nThe prime factorization of the conductor is d = {d} = {' * '.join(map(str, prime_factors))}.")
    
    if any(exp > 1 for exp in d_factors.values()):
        print("The conductor is not square-free. This script is designed for square-free conductors.")
        return

    # Step 2: Set up the calculation using Mobius inversion.
    print(f"We need to find the number of characters with order exactly {k}.")
    print("The formula is: Sum_{n|k} mu(k/n) * g(n)")
    print("where g(n) is the number of tuples of primitive characters (chi_i) where each ord(chi_i) divides n.")
    
    p_minus_1_list = [p - 1 for p in prime_factors]
    
    def g(n, p_minus_1s):
        """
        Calculates the number of character tuples with order dividing n.
        """
        if n == 0: return 0
        count = 1
        for p_m_1 in p_minus_1s:
            # Number of non-trivial characters mod p with order dividing n
            num_chars = math.gcd(n, p_m_1) - 1
            count *= num_chars
        return count

    # Step 3: Perform the calculation.
    divisors_of_k = get_divisors(k)
    total = 0
    
    print("\nCalculation steps:")
    value_calc_parts = []
    
    for n in divisors_of_k:
        mu_val = mobius(k // n)
        g_val = g(n, p_minus_1_list)
        term = mu_val * g_val
        total += term
        
        op = "+" if term >= 0 else "-"
        # Build the final equation string part by part
        if term >= 0:
            value_calc_parts.append(f"+ {term}")
        else:
            value_calc_parts.append(f"- {abs(term)}")
        
        print(f"  For n={n}: mu({k//n})*g({n}) = {mu_val} * {g_val} = {term}")

    # Step 4: Display the final result and equation.
    # We clean up the leading '+' sign for the equation string.
    final_equation = " ".join(value_calc_parts).lstrip("+ ")

    print("\nThe final calculation is:")
    print(f"Result = {final_equation}")
    print(f"Result = {total}")

if __name__ == '__main__':
    main()