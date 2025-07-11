import math

def phi(n):
    """
    Calculates Euler's totient function.
    """
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
    """
    Returns the prime factorization of a number as a dictionary.
    """
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

def get_divisors(n):
    """
    Returns a sorted list of divisors of a number.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

# Set the conductor N and the desired order
N = 36036
ORDER = 6

print(f"To find the number of primitive Dirichlet characters of conductor N = {N} and order {ORDER}, we first factor N into prime powers.")
prime_factors = get_prime_factorization(N)
# Sort factors for consistent output order, e.g., by p
sorted_factors = sorted(prime_factors.items())
factor_str = ' * '.join([f'{p}^{a}' for p, a in sorted_factors])
print(f"The factorization is N = {factor_str}.")
print("A character is primitive modulo N if and only if its component characters are primitive modulo each prime power factor.")
print(f"We count the number of valid primitive characters for each factor whose order divides {ORDER}.\n")

total_count = 1
counts_per_factor = []
factor_strings = []

for p, a in sorted_factors:
    pa = p**a
    
    # Determine possible orders of primitive characters for p^a
    if p == 2:
        if a == 2: # Conductor is 4
            # Primitive character mod 4 has order 2.
            primitive_orders = [2]
        else: # Other powers of 2 (not relevant for N=36036)
            primitive_orders = []
    else: # p is an odd prime
        phi_pa = phi(pa)
        divisors_phi_pa = get_divisors(phi_pa)
        if a > 1:
            phi_prev = phi(p**(a-1))
        else:
            phi_prev = phi(1)
        
        divisors_phi_prev = get_divisors(phi_prev)
        # A character is primitive if its order does not divide phi of the next lower power of p.
        primitive_orders = [d for d in divisors_phi_pa if d not in divisors_phi_prev]

    # Filter primitive orders to those dividing the required ORDER
    valid_orders = [d for d in primitive_orders if ORDER % d == 0]
    
    count_for_factor = 0
    for d in valid_orders:
        count_for_factor += phi(d)
        
    counts_per_factor.append(count_for_factor)
    
    # Prepare explanation string for this factor
    calc_str = " + ".join([f"phi({d})" for d in valid_orders]) if valid_orders else "0"
    val_str = " + ".join([str(phi(d)) for d in valid_orders]) if valid_orders else "0"
    
    factor_strings.append(f"For the factor {pa}: The number of primitive characters with order dividing {ORDER} is {calc_str} = {val_str} = {count_for_factor}.")

for s in factor_strings:
    print(s)

print(f"\nThe lcm of the orders is guaranteed to be {ORDER} because the character choices for factors 4 and 9 ensure divisibility by 2 and 3 respectively.")
print("The total number of such characters is the product of the counts for each factor.")
final_calc_str = " * ".join(map(str, counts_per_factor))
print(f"Total number = {final_calc_str} = {math.prod(counts_per_factor)}")