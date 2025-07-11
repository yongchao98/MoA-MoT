import math

def phi(n):
    """
    Computes Euler's totient function phi(n).
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
    Calculates the number of primitive Dirichlet characters of conductor N=36036 and order 6.
    """
    N = 36036
    print(f"We want to find the number of primitive Dirichlet characters of conductor N = {N} and order 6.")
    
    print("\nStep 1: Factor the conductor N into its prime power components.")
    factors_str = "2^2 * 3^2 * 7 * 11 * 13"
    factors_vals = "4 * 9 * 7 * 11 * 13"
    print(f"N = 36036 = {factors_str} = {factors_vals}\n")

    print("Step 2: For each prime power factor q, find the number of primitive characters modulo q")
    print("whose order divides 6.\n")

    # For q = 4
    # The group of characters mod 4 is cyclic of order phi(4)=2. The single primitive character has order 2.
    choices_4 = phi(2)
    print(f"For q = 4: The group of characters is C_2. Primitive characters have order 2.")
    print(f"  Number of choices with order dividing 6 is phi(2) = {choices_4}.\n")

    # For q = 9
    # The group of characters mod 9 is cyclic of order phi(9)=6. Primitive characters have orders 3 and 6.
    choices_9 = phi(3) + phi(6)
    print(f"For q = 9: The group of characters is C_6. Primitive characters have orders 3 and 6.")
    print(f"  Number of choices with order dividing 6 is phi(3) + phi(6) = {phi(3)} + {phi(6)} = {choices_9}.\n")

    # For q = 7
    # The group of characters mod 7 is cyclic of order phi(7)=6. Primitive characters have orders 2, 3, 6.
    choices_7 = phi(2) + phi(3) + phi(6)
    print(f"For q = 7: The group of characters is C_6. Primitive characters have orders 2, 3, 6.")
    print(f"  Number of choices with order dividing 6 is phi(2) + phi(3) + phi(6) = {phi(2)} + {phi(3)} + {phi(6)} = {choices_7}.\n")
    
    # For q = 11
    # The group of characters mod 11 is C_10. We need primitive characters with order dividing 6. Only order 2 works.
    choices_11 = phi(2)
    print(f"For q = 11: The group of characters is C_10. Primitive characters whose orders divide 6 can only have order 2.")
    print(f"  Number of choices is phi(2) = {choices_11}.\n")

    # For q = 13
    # The group of characters mod 13 is C_12. Primitive characters with order dividing 6 are 2, 3, 6.
    choices_13 = phi(2) + phi(3) + phi(6)
    print(f"For q = 13: The group of characters is C_12. Primitive characters whose orders divide 6 can have order 2, 3, or 6.")
    print(f"  Number of choices is phi(2) + phi(3) + phi(6) = {phi(2)} + {phi(3)} + {phi(6)} = {choices_13}.\n")
    
    print("Step 3: Combine the results.")
    print("The order of a character is the lcm of the orders of its components. Due to the factors 4 and 9,")
    print("the final order must be a multiple of lcm(2, 3) = 6. As all component orders divide 6,")
    print("the final order is guaranteed to be exactly 6.\n")
    
    print("The total number of characters is the product of the number of choices for each factor:")
    total = choices_4 * choices_9 * choices_7 * choices_11 * choices_13
    
    # Final equation printout
    print(f"{choices_4} (from q=4) * {choices_9} (from q=9) * {choices_7} (from q=7) * {choices_11} (from q=11) * {choices_13} (from q=13) = {total}")

if __name__ == '__main__':
    solve()
    print("\n<<<100>>>")
