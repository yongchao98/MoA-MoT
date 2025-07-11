import math

def phi(n):
    """
    Calculates Euler's totient function phi(n).
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

def solve_and_print():
    """
    Solves the problem and prints the detailed steps and final answer.
    """
    N = 36036
    order = 6
    
    print(f"We want to find the number of primitive Dirichlet characters of conductor N = {N} and order {order}.")
    print("-" * 70)
    
    # Step 1: Prime factorization of N
    factors = [4, 9, 7, 11, 13]
    print(f"Step 1: The prime factorization of N is 2^2 * 3^2 * 7 * 11 * 13.")
    print("A primitive character chi of conductor N can be written as a product:")
    print("chi = chi_4 * chi_9 * chi_7 * chi_11 * chi_13, where each chi_m is a primitive character modulo m.")
    print(f"The order of chi, lcm(ord(chi_4), ..., ord(chi_13)), must be {order}.")
    print("-" * 70)
    
    print(f"Step 2: For each modulus, find the number of primitive characters whose order divides {order}.")

    # Modulo 4
    count_4 = phi(2)
    print("\nFor modulus 4 (2^2):")
    print("  - Primitive characters mod 4 must have order 2.")
    print("  - The order 2 divides the required order 6.")
    print(f"  - Number of such characters = phi(2) = {count_4}")
    
    # Modulo 9
    count_9 = phi(3) + phi(6)
    print("\nFor modulus 9 (3^2):")
    print("  - Primitive characters mod 9 must have order 3 or 6.")
    print("  - Both orders 3 and 6 divide the required order 6.")
    print(f"  - Number of such characters = phi(3) + phi(6) = {phi(3)} + {phi(6)} = {count_9}")

    # Modulo 7
    count_7 = phi(2) + phi(3) + phi(6)
    print("\nFor modulus 7:")
    print("  - Primitive characters mod 7 have orders 2, 3, or 6.")
    print("  - All these orders divide the required order 6.")
    print(f"  - Number of such characters = phi(2) + phi(3) + phi(6) = {phi(2)} + {phi(3)} + {phi(6)} = {count_7}")
    
    # Modulo 11
    count_11 = phi(2)
    print("\nFor modulus 11:")
    print("  - Primitive characters mod 11 have orders 2, 5, or 10.")
    print("  - Of these, only order 2 divides the required order 6.")
    print(f"  - Number of such characters = phi(2) = {phi(2)}")

    # Modulo 13
    count_13 = phi(2) + phi(3) + phi(6)
    print("\nFor modulus 13:")
    print("  - Primitive characters mod 13 have orders 2, 3, 4, 6, or 12.")
    print("  - Of these, orders 2, 3, and 6 divide the required order 6.")
    print(f"  - Number of such characters = phi(2) + phi(3) + phi(6) = {phi(2)} + {phi(3)} + {phi(6)} = {count_13}")
    print("-" * 70)

    # Step 3: Calculate the total number and verify the order
    total_count = count_4 * count_9 * count_7 * count_11 * count_13
    print("Step 3: The total number is the product of the counts for each modulus.")
    print("We must confirm the resulting character's order is exactly 6.")
    print("  - The order is a divisor of 6 because all component orders divide 6.")
    print("  - The order is a multiple of 2 because the component from mod 4 has order 2.")
    print("  - The order is a multiple of 3 because the component from mod 9 has order 3 or 6.")
    print("An integer that is a multiple of 6 and a divisor of 6 must be exactly 6.")
    print("-" * 70)
    
    print("Final Calculation:")
    print("The number of primitive characters of conductor 36036 and order 6 is:")
    print(f"{count_4} (from mod 4) * {count_9} (from mod 9) * {count_7} (from mod 7) * {count_11} (from mod 11) * {count_13} (from mod 13) = {total_count}")

solve_and_print()