import math

def phi(n):
    """
    Calculates Euler's totient function phi(n), which counts the number of
    positive integers up to a given integer n that are relatively prime to n.
    This is used to find the number of characters of a specific order.
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

def solve():
    """
    Calculates the number of primitive Dirichlet characters of a given
    conductor and order by multiplying the number of choices for each
    component character derived from the prime factorization of the conductor.
    """
    N = 36036
    print(f"We want to find the number of primitive Dirichlet characters of conductor N = {N} and order 6.")
    
    print("\nStep 1: The prime factorization of N is 2^2 * 3^2 * 7 * 11 * 13.")
    
    # Number of choices for the component character modulo 4
    # There is 1 primitive character mod 4, and its order is 2. Since 2 divides 6, there is 1 choice.
    choices_4 = 1
    print(f"\nFor modulus 4 (2^2), there is {choices_4} choice.")

    # Number of choices for the component character modulo 9
    # Primitive characters have orders 3 or 6. Both divide 6.
    choices_9 = phi(3) + phi(6)
    print(f"For modulus 9 (3^2), the number of choices is phi(3) + phi(6) = {phi(3)} + {phi(6)} = {choices_9}.")

    # Number of choices for the component character modulo 7
    # Primitive characters have orders 2, 3, or 6. All divide 6.
    choices_7 = phi(2) + phi(3) + phi(6)
    print(f"For modulus 7, the number of choices is phi(2) + phi(3) + phi(6) = {phi(2)} + {phi(3)} + {phi(6)} = {choices_7}.")
    
    # Number of choices for the component character modulo 11
    # Primitive characters have orders 2, 5, or 10. Only order 2 divides 6.
    choices_11 = phi(2)
    print(f"For modulus 11, the number of choices is phi(2) = {choices_11}.")

    # Number of choices for the component character modulo 13
    # Primitive characters have orders 2, 3, 4, 6, 12. Orders 2, 3, 6 divide 6.
    choices_13 = phi(2) + phi(3) + phi(6)
    print(f"For modulus 13, the number of choices is phi(2) + phi(3) + phi(6) = {phi(2)} + {phi(3)} + {phi(6)} = {choices_13}.")
    
    # The total number is the product of the number of choices for each component.
    total = choices_4 * choices_9 * choices_7 * choices_11 * choices_13
    
    print("\nStep 2: The total number is the product of these choices.")
    print(f"Total number = {choices_4} * {choices_9} * {choices_7} * {choices_11} * {choices_13}")
    print(f"Total number = {total}")

solve()