import math
from sympy import factorint, totient

def num_primitive_chars_prime_power(p, a, k):
    """
    Calculates the number of primitive Dirichlet characters of conductor p^a and order k.
    """
    if k == 0:
        return 0
    # The trivial character (order 1) is only primitive for conductor 1.
    if k == 1:
        return 0

    # Case p=2
    if p == 2:
        if a == 1: # Conductor 2 does not have primitive characters. The only char is trivial mod 1.
            return 0
        if a == 2: # Conductor 4
            return 1 if k == 2 else 0
        if a == 3: # Conductor 8
            # The group of characters is C2 x C2. Both primitive characters have order 2.
            return 2 if k == 2 else 0
        # For a > 3, not needed for this problem, but formula gets more complex.

    # Case p is an odd prime
    if p % 2 != 0:
        if a == 1:
            # For a prime conductor p, all non-trivial characters are primitive.
            # The group of characters is cyclic of order p-1.
            if k % (p - 1) == 0:
                return totient(k)
            else:
                return 0
        else: # a > 1
            # For p^a, primitive characters are those not induced from p^(a-1).
            # The group of chars mod p^a is cyclic of order phi(p^a).
            # A char of order k is primitive iff k | phi(p^a) but k does not divide phi(p^(a-1)).
            phi_pa = totient(p**a)
            phi_pa_minus_1 = totient(p**(a-1))
            if k % phi_pa == 0 and k % phi_pa_minus_1 != 0:
                return totient(k)
            else:
                return 0
    return 0

def solve():
    """
    Solves the problem of finding the number of primitive Dirichlet characters
    of conductor N=36036 and order 6.
    """
    N = 36036
    factors = factorint(N)

    print(f"The conductor is N = 36036.")
    factor_str = " * ".join([f"{p}^{a}" for p, a in factors.items()])
    print(f"Step 1: Factor the conductor N = {factor_str}")
    
    print("\nStep 2: A character is primitive of conductor N if it is a product of primitive characters for each prime power factor.")
    print("The order of the character is the LCM of the orders of its components.")
    print("We need the LCM of the orders to be 6.")

    print("\nStep 3: For each factor, find the number of primitive characters with order dividing 6.")
    
    conductors = {p**a: (p, a) for p, a in factors.items()}
    choices = {}

    # Conductor 4 (2^2)
    p, a = conductors[4]
    n4 = num_primitive_chars_prime_power(p, a, 2)
    choices[4] = n4
    print(f"\nFor conductor 4 = 2^2:")
    print(f"  Primitive characters must have order 2. Number of such characters: {n4}.")
    print(f"  So, for the component modulo 4, there is {n4} choice, and its order must be 2.")
    
    # Conductor 9 (3^2)
    p, a = conductors[9]
    n9_3 = num_primitive_chars_prime_power(p, a, 3)
    n9_6 = num_primitive_chars_prime_power(p, a, 6)
    choices[9] = n9_3 + n9_6
    print(f"\nFor conductor 9 = 3^2:")
    print(f"  Primitive characters can have order 3 (number: {n9_3}) or 6 (number: {n9_6}).")
    print(f"  Total choices for component modulo 9, with order dividing 6: {n9_3} + {n9_6} = {choices[9]}.")
    print(f"  Note: The order is always a multiple of 3.")

    # Conductor 7
    p, a = conductors[7]
    n7_2 = num_primitive_chars_prime_power(p, a, 2)
    n7_3 = num_primitive_chars_prime_power(p, a, 3)
    n7_6 = num_primitive_chars_prime_power(p, a, 6)
    choices[7] = n7_2 + n7_3 + n7_6
    print(f"\nFor conductor 7:")
    print(f"  Primitive characters can have order 2 (num: {n7_2}), 3 (num: {n7_3}), or 6 (num: {n7_6}).")
    print(f"  Total choices for component modulo 7, with order dividing 6: {n7_2} + {n7_3} + {n7_6} = {choices[7]}.")

    # Conductor 11
    p, a = conductors[11]
    # Order must divide gcd(6, phi(11)=10) = 2.
    n11_2 = num_primitive_chars_prime_power(p, a, 2)
    choices[11] = n11_2
    print(f"\nFor conductor 11:")
    print(f"  The order must divide gcd(6, 10)=2. So, it must be 2.")
    print(f"  Number of primitive characters of order 2: {n11_2}.")
    print(f"  Total choices for component modulo 11: {choices[11]}. Its order must be 2.")
    
    # Conductor 13
    p, a = conductors[13]
    # Order must divide gcd(6, phi(13)=12) = 6.
    n13_2 = num_primitive_chars_prime_power(p, a, 2)
    n13_3 = num_primitive_chars_prime_power(p, a, 3)
    n13_6 = num_primitive_chars_prime_power(p, a, 6)
    choices[13] = n13_2 + n13_3 + n13_6
    print(f"\nFor conductor 13:")
    print(f"  The order must divide gcd(6, 12)=6. It can be 2, 3, or 6.")
    print(f"  Number of characters: order 2 ({n13_2}), order 3 ({n13_3}), order 6 ({n13_6}).")
    print(f"  Total choices for component modulo 13: {n13_2} + {n13_3} + {n13_6} = {choices[13]}.")

    print("\nStep 4: Combine the results.")
    print("The final character's order is the LCM of component orders.")
    print("Let the orders be (k_4, k_9, k_7, k_11, k_13).")
    print("- At least one order must be a multiple of 2 (guaranteed by conductors 4 and 11).")
    print("- At least one order must be a multiple of 3 (guaranteed by conductor 9).")
    print("- All orders are divisors of 6.")
    print("Therefore, the LCM is guaranteed to be 6 for any combination of the chosen characters.")
    
    total = 1
    for c in choices.values():
        total *= c

    print("\nThe total number is the product of the number of choices for each conductor:")
    print(f"Total = {choices[4]} (for 4) * {choices[9]} (for 9) * {choices[7]} (for 7) * {choices[11]} (for 11) * {choices[13]} (for 13)")
    print(f"Total = {choices[4]} * {choices[9]} * {choices[7]} * {choices[11]} * {choices[13]} = {total}")

    return total

if __name__ == '__main__':
    result = solve()
    # The final answer in the requested format will be printed after the explanation
    # For now, just ensuring the logic holds. This final part won't be in the response.
    # print(f"\n<<< {result} >>>")