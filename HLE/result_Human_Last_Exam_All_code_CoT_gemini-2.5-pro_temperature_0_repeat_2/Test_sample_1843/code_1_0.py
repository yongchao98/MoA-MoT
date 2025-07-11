import math

def phi(n):
    """
    Calculates Euler's totient function, phi(n), which counts the number of
    positive integers up to a given integer n that are relatively prime to n.
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
    return result

def solve():
    """
    Calculates the number of primitive Dirichlet characters of conductor 36036 and order 6.
    """
    N = 36036
    print(f"We want to find the number of primitive Dirichlet characters of conductor N = {N} and order 6.")
    print("-" * 80)
    
    print("Step 1: Factorize the conductor N.")
    # N = 2^2 * 3^2 * 7 * 11 * 13
    print(f"The prime factorization of N = {N} is 2^2 * 3^2 * 7 * 11 * 13.")
    print("A primitive character chi mod N can be uniquely decomposed as a product:")
    print("chi = chi_4 * chi_9 * chi_7 * chi_11 * chi_13")
    print("where each chi_{p^a} is a primitive character modulo the prime power factor p^a.")
    print("-" * 80)

    print("Step 2: Analyze the order condition.")
    print("The order of chi is lcm(ord(chi_4), ord(chi_9), ord(chi_7), ord(chi_11), ord(chi_13)).")
    print("For the order to be 6, the order of each component chi_{p^a} must divide 6.")
    print("-" * 80)

    print("Step 3: Count the number of valid primitive components for each prime power factor.")
    
    # Conductor 4 = 2^2
    print("\nFor conductor 4 (2^2):")
    num_prim_chars_4 = phi(4) - phi(2)
    print(f"The number of primitive characters is phi(4) - phi(2) = {phi(4)} - {phi(2)} = {num_prim_chars_4}.")
    print("This character has order 2. Since 2 divides 6, it is a valid component.")
    c4 = 1
    print(f"Number of choices for chi_4: {c4}.")

    # Conductor 9 = 3^2
    print("\nFor conductor 9 (3^2):")
    num_prim_chars_9 = phi(9) - phi(3)
    print(f"The number of primitive characters is phi(9) - phi(3) = {phi(9)} - {phi(3)} = {num_prim_chars_9}.")
    print("Their orders are 3 and 6. All these orders divide 6.")
    c9 = num_prim_chars_9
    print(f"Number of choices for chi_9: {c9}.")

    # Conductor 7
    print("\nFor conductor 7:")
    num_prim_chars_7 = phi(7) - 1
    print(f"The primitive characters are the non-principal ones, numbering phi(7) - 1 = {num_prim_chars_7}.")
    print("Their orders are 2, 3, and 6. All these orders divide 6.")
    c7 = num_prim_chars_7
    print(f"Number of choices for chi_7: {c7}.")

    # Conductor 11
    print("\nFor conductor 11:")
    print("The group of characters mod 11 is cyclic of order phi(11) = 10.")
    print("Primitive characters have orders that divide 10 (except 1), i.e., 2, 5, 10.")
    print("For the order to divide 6, it must be 2.")
    c11 = phi(2)
    print(f"The number of characters of order 2 is phi(2) = {c11}.")
    print(f"Number of choices for chi_11: {c11}.")

    # Conductor 13
    print("\nFor conductor 13:")
    print("The group of characters mod 13 is cyclic of order phi(13) = 12.")
    print("Primitive characters have orders that divide 12 (except 1).")
    print("For the order to divide 6, it can be 2, 3, or 6.")
    c13_ord2 = phi(2)
    c13_ord3 = phi(3)
    c13_ord6 = phi(6)
    c13 = c13_ord2 + c13_ord3 + c13_ord6
    print(f"Number of choices for chi_13 = (chars of order 2) + (chars of order 3) + (chars of order 6)")
    print(f"  = phi(2) + phi(3) + phi(6) = {c13_ord2} + {c13_ord3} + {c13_ord6} = {c13}.")
    
    print("-" * 80)
    print("Step 4: Verify the lcm condition and calculate the total number.")
    print("For any combination of the chosen components:")
    print("- The order of chi_4 is 2, so the lcm of all orders is a multiple of 2.")
    print("- The order of chi_9 is 3 or 6, so the lcm is a multiple of 3.")
    print("Since the lcm is a multiple of 2 and 3, it is a multiple of 6.")
    print("As all component orders divide 6, the lcm must be exactly 6.")
    print("\nThe total number is the product of the number of choices for each component.")
    
    total = c4 * c9 * c7 * c11 * c13
    print("\nFinal Calculation:")
    print(f"{c4} (from 4) * {c9} (from 9) * {c7} (from 7) * {c11} (from 11) * {c13} (from 13) = {total}")

solve()
<<<100>>>