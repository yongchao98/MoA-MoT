import math
from collections import defaultdict

def prime_factorize(n):
    """Computes the prime factorization of a number n."""
    factors = defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return dict(factors)

def phi(n):
    """Computes Euler's totient function."""
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

def solve():
    """
    Calculates the number of primitive Dirichlet characters for a given conductor N and order.
    """
    N = 36036
    ORDER = 6
    print(f"We want to find the number of primitive Dirichlet characters of conductor N = {N} and order {ORDER}.")
    print("-" * 80)

    # Step 1: Prime factorization of N
    factors = prime_factorize(N)
    factor_str = " * ".join([f"{p}^{a}" if a > 1 else str(p) for p, a in sorted(factors.items())])
    print("Step 1: Prime factorize the conductor N.")
    print(f"N = {N} = {factor_str}")
    print("-" * 80)

    # Step 2: Character decomposition
    print("Step 2: Decompose the character based on the factorization of N.")
    print("A character chi mod N can be uniquely written as a product of characters: chi = chi_4 * chi_9 * chi_7 * chi_11 * chi_13.")
    print("chi is primitive with conductor N if and only if each component chi_p^a is primitive.")
    print(f"The order of chi, lcm(ord(chi_4), ord(chi_9), ...), must be {ORDER}.")
    print("-" * 80)
    
    # Step 3: Analyze each component
    print("Step 3: Count the valid primitive characters for each prime power factor.")
    
    num_choices = {}
    
    # Iterate through prime factors p^a of N
    for p, a in sorted(factors.items()):
        m = p**a
        
        # Determine the structure of the character group
        group_order = phi(m)
        
        primitive_orders_info = []
        count = 0

        # Check all possible orders 'd' that a character in this group can have.
        for d in range(1, group_order + 1):
            if group_order % d == 0: # 'd' must divide the group order
                # Condition 1: The order 'd' of the component must divide the target order
                if ORDER % d == 0:
                    is_primitive = False
                    # Condition 2: The character must be primitive for its modulus
                    if a == 1: # For p^1, primitive means not the principal character (order > 1)
                        if d > 1:
                            is_primitive = True
                    elif p == 2 and a == 2: # For 2^2=4, primitive means not principal (order > 1)
                        if d > 1:
                            is_primitive = True
                    elif p%2 == 1 and a > 1: # For p^a, p odd, primitive means order 'd' does not divide phi(p^(a-1))
                        if phi(p**(a-1)) % d != 0:
                           is_primitive = True

                    if is_primitive:
                        # Number of elements of order d in a cyclic group is phi(d)
                        num_d = phi(d)
                        primitive_orders_info.append(f"{num_d} of order {d}")
                        count += num_d
        
        num_choices[m] = count
        
        print(f"For modulus {m} = {p}^{a}:")
        print(f"  - The character group is cyclic of order phi({m}) = {group_order}.")
        print(f"  - We need primitive characters whose order divides {ORDER}.")
        if primitive_orders_info:
            print(f"  - Available primitive characters: {' and '.join(primitive_orders_info)}.")
            print(f"  - Total choices for mod {m}: {count}.")
        else:
            print(f"  - No such characters exist.")

    print("-" * 80)
    
    # Step 4: Final calculation
    print("Step 4: Combine the results.")
    print("The order of the final character is the lcm of the component orders.")
    print(f"For the final order to be {ORDER}, the lcm must be a multiple of 2 and 3, and all component orders must divide {ORDER}.")
    print(" - The character mod 4 must have order 2. This guarantees the lcm is a multiple of 2.")
    print(" - The character mod 9 must have order 3 or 6. This guarantees the lcm is a multiple of 3.")
    print(f" - All allowed component character orders found above are divisors of {ORDER}.")
    print(f"Therefore, any combination of these selected components will have an order of exactly {ORDER}.")
    print("The total number of such characters is the product of the number of choices for each component.")
    
    # Final equation
    equation_parts = []
    final_product = 1
    for m in sorted(num_choices.keys()):
        equation_parts.append(str(num_choices[m]))
        final_product *= num_choices[m]

    print()
    print("Final Calculation:")
    print(f"Number of characters = { ' * '.join(equation_parts) } = {final_product}")

solve()
<<<100>>>