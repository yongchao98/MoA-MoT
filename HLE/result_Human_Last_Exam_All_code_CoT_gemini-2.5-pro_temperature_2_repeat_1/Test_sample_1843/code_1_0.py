import math

def phi(n):
    """Euler's totient function."""
    if n == 1:
        return 1
    result = n
    p = 2
    temp_n = n
    while p * p <= temp_n:
        if temp_n % p == 0:
            while temp_n % p == 0:
                temp_n //= p
            result -= result // p
        p += 1
    if temp_n > 1:
        result -= result // temp_n
    return int(result)

def get_prime_factorization(n):
    """Returns a dictionary of prime factors and their powers."""
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def gcd(a, b):
    """Greatest common divisor."""
    return math.gcd(a, b)

def count_primitive_chars_by_order_divisor(p, k, D):
    """
    Counts the number of primitive Dirichlet characters modulo p^k
    whose order divides D. Works for p odd, or p^k=4.
    """
    val_pk = phi(p**k)
    # phi(p^0) represents the character group mod 1, which has 1 element
    val_pk_minus_1 = phi(p**(k - 1)) if k > 1 else 1
    return gcd(val_pk, D) - gcd(val_pk_minus_1, D)

def solve():
    """
    Main function to solve the problem.
    """
    N = 36036
    ORDER = 6

    print(f"To find the number of primitive Dirichlet characters of conductor N = {N} and order {ORDER}, we proceed as follows:")
    
    # Step 1: Factor N
    factors = get_prime_factorization(N)
    factor_str = " * ".join([f"{p}^{k}" if k > 1 else str(p) for p, k in sorted(factors.items())])
    print(f"\n1. The prime factorization of N is: {N} = {factor_str}.")

    # Step 2: Count choices for each factor
    print("\n2. We count the number of primitive characters for each prime power factor whose order divides 6:")
    counts = {}
    for p, k in sorted(factors.items()):
        count = count_primitive_chars_by_order_divisor(p, k, ORDER)
        counts[(p,k)] = count
        modulus_str = f"{p}^{k}" if k > 1 else str(p)
        print(f"   - For modulus {modulus_str}: There are {count} such characters.")
        
    # Step 3: Check LCM condition
    print("\n3. We verify that the final order is exactly 6:")
    print("   - For modulus 4, any primitive character must have order 2. This makes the lcm divisible by 2.")
    print("   - For modulus 9, any primitive character must have order 3 or 6. This makes the lcm divisible by 3.")
    print("   Since the lcm of the orders must be divisible by 2 and 3, it must be divisible by 6.")
    print("   As all component orders divide 6, their lcm must be exactly 6.")
    
    # Step 4: Calculate total number
    total_count = 1
    count_list = []
    for p_k in sorted(factors.items()):
        p, k = p_k
        count = counts[(p, k)]
        total_count *= count
        count_list.append(str(count))
            
    print("\n4. The total number of characters is the product of the number of choices for each factor.")
    equation = " * ".join(count_list)
    print(f"   Number = {equation} = {total_count}")

solve()