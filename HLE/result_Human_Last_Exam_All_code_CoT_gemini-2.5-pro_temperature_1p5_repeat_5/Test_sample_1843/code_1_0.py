import math

def phi(n):
    """
    Computes Euler's totient function phi(n).
    """
    if n == 1:
        return 1
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

def count_primitive_chars_by_order(modulus, required_order_lcm):
    """
    Counts the number of primitive characters for a given modulus whose orders
    divide required_order_lcm.
    """
    # Orders of primitive characters for the given moduli
    # whose orders divide 6
    
    count = 0
    if modulus == 4:
        # One primitive char of order 2. 2 divides 6.
        count = 1
    elif modulus == 9:
        # Primitive chars have orders 3 and 6. Both divide 6.
        # Number of primitive chars of order k is phi(k).
        count = phi(3) + phi(6)
    elif modulus in [7, 13]:
        # For a prime p, primitive chars are non-principal ones.
        # Their orders are divisors of p-1 (excluding 1).
        # We need orders that also divide 6.
        if modulus == 7: # p-1 = 6
            # orders are 2, 3, 6. All divide 6.
            count = phi(2) + phi(3) + phi(6)
        if modulus == 13: # p-1 = 12
            # orders are 2, 3, 4, 6, 12. Those dividing 6 are 2, 3, 6.
            count = phi(2) + phi(3) + phi(6)
    elif modulus == 11:
        # p-1 = 10. orders are 2, 5, 10. Only 2 divides 6.
        count = phi(2)
        
    return count
    
def solve():
    """
    Solves the problem of finding the number of primitive Dirichlet characters
    of conductor N=36036 and order 6.
    """
    N = 36036
    order = 6
    
    print(f"To find the number of primitive Dirichlet characters of conductor N = {N} and order {order}, we first factorize N.")
    print(f"N = 36036 = 2^2 * 3^2 * 7 * 11 * 13.")
    
    moduli = [4, 9, 7, 11, 13]
    
    print("\nA character is primitive with conductor N if and only if it is a product of primitive characters for each prime power factor.")
    print(f"The order of the character is the lcm of the orders of its components. We need this lcm to be {order}.")
    print("\nThis means the order of each component character must divide 6. We count the number of choices for each component:")

    counts = []
    for m in moduli:
        c = count_primitive_chars_by_order(m, order)
        counts.append(c)
        print(f"- For modulus {m}, there are {c} primitive characters whose order divides 6.")

    total_count = math.prod(counts)

    print("\nThe total number of such combinations of characters is the product of the individual counts.")
    
    # Building the equation string
    equation_parts = [str(c) for c in counts]
    equation_str = " * ".join(equation_parts)
    
    print(f"\nFinal calculation: {equation_str} = {total_count}")

    print("\nWe must verify the order is exactly 6. The order of the component for modulus 4 is always 2, and the order for modulus 9 is a multiple of 3.")
    print("Therefore, the lcm of the orders is always a multiple of lcm(2, 3) = 6.")
    print("Since we already ensured the component orders divide 6, the final order must be exactly 6.")
    print(f"So, the total number of such characters is {total_count}.")


solve()
<<<100>>>