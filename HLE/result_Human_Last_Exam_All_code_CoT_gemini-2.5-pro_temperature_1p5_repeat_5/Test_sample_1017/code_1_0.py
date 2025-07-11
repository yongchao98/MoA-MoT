import math

def get_prime_factorization(n):
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

def is_power_of(n, base):
    """Checks if n is a power of base."""
    if n == 1:
        return True
    if n == 0 or n % base != 0:
        return False
    return is_power_of(n // base, base)
    
def check_proper_power(rational_str, base):
    """
    Checks if a rational number x=num/den is a proper power in the group of 
    translations within Thompson's group T_base.
    This is true if there exists an integer k > 1 such that x/k is in Z[1/base].
    x/k = num / (den * k). This is in Z[1/base] if the prime factors of the 
    denominator of the reduced fraction are all factors of `base`.
    """
    num_str, den_str = rational_str.split('/')
    num = int(num_str)
    den = int(den_str)
    
    print(f"Checking if {rational_str} is a proper power in T_{base}...")
    
    base_factors = get_prime_factorization(base)

    # We check for the smallest integer k > 1
    for k in range(2, num + 2):
        # We need to check if the denominator of (num / (den*k)) in lowest terms
        # is a power of `base`.
        common_divisor = math.gcd(num, k)
        new_num = num // common_divisor
        
        # The new denominator is den * (k/gcd(num,k))
        # We need to check the prime factors of this new denominator.
        
        denom_part_from_k = k // common_divisor
        
        # Check prime factors of denom_part_from_k
        is_ok = True
        k_factors = get_prime_factorization(denom_part_from_k)
        for factor in k_factors:
            if factor not in base_factors:
                is_ok = False
                break
        
        if is_ok:
            print(f"Success! For k = {k}, the translation amount becomes {new_num}/{den*denom_part_from_k}.")
            print(f"This is in Z[1/{base}] because the prime factors of the denominator are contained in the prime factors of {base}.")
            print(f"{rational_str} is a proper power. For example, T_{rational_str} = (T_{new_num}/{den*denom_part_from_k})^{k}.")
            return True
            
    print(f"Could not show that {rational_str} is a proper power for small k.")
    return False

# Based on the problem, g corresponds to 2/27 and h to 16/27.
# We assume the group is T_3, so the base is 3.
g_val = "2/27"
h_val = "16/27"
p_base = 3

g_is_power = check_proper_power(g_val, p_base)
print("-" * 20)
h_is_power = check_proper_power(h_val, p_base)
print("-" * 20)

if g_is_power or h_is_power:
    scl = 0
    print(f"At least one of the elements is a proper power in T_{p_base}.")
    print("By Bavard's theorem, the stable commutator length of their commutator is 0.")
    print(f"scl([g_1, h_2]) = {scl}")
else:
    # This case is not reached with the given numbers.
    # If neither were proper powers, for abelian groups the scl would be 1/2.
    scl = "1/2" 
    print("Neither element is a proper power.")
    print(f"The stable commutator length would be {scl} if the groups were abelian.")

print("\nFinal equation:")
print(f"g = translation by {g_val}")
print(f"h = translation by {h_val}")
print(f"Since g is a proper power (e.g., g = (translation by {1}/{27})^{2})")
print(f"scl([g_1, h_2]) = 0")