import math

def phi(n):
    """Computes Euler's totient function."""
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

def solve_character_count():
    """
    Calculates the number of primitive Dirichlet characters of conductor N=36036 and order 6.
    """
    N = 36036
    print(f"We want to find the number of primitive Dirichlet characters of conductor N = {N} and order 6.")
    
    # Step 1: Prime factorization of N
    print("\nStep 1: Prime factorization of N")
    factorization = "2^2 * 3^2 * 7 * 11 * 13"
    print(f"N = 36036 = {factorization}")
    
    print("\nA character chi mod N is primitive if it's a product of primitive characters for each prime power factor.")
    print("The order of chi is the lcm of the orders of its component characters.")
    print("The order of each component character must divide 6.\n")
    
    # Step 2: Analyze each prime power factor
    
    # For N = 4 = 2^2
    print("--- Component for conductor 4 (2^2) ---")
    group_structure_4 = "C_2"
    num_prim_4 = 1 # The non-trivial character mod 4
    order_4 = 2
    print(f"The group of characters mod 4 is isomorphic to {group_structure_4}. It has {phi(4)} members.")
    print("A character is primitive mod 4 if its conductor is 4. There is 1 such character.")
    print(f"Its order is {order_4}, which divides 6. So it is a valid choice.")
    choices_4 = num_prim_4
    print(f"Number of choices for conductor 4: {choices_4}")

    # For N = 9 = 3^2
    print("\n--- Component for conductor 9 (3^2) ---")
    group_structure_9 = "C_6"
    num_prim_9 = phi(9) - phi(3)
    print(f"The group of characters mod 9 is isomorphic to {group_structure_9}. It has {phi(9)} members.")
    print(f"The number of primitive characters mod 9 is phi(9) - phi(3) = {phi(9)} - {phi(3)} = {num_prim_9}.")
    # Orders of primitive characters mod 9 are 3 and 6
    print("The orders of these primitive characters are 3 (two of them) and 6 (two of them).")
    print("All these orders (3, 6) divide 6. So all 4 primitive characters are valid choices.")
    choices_9 = num_prim_9
    print(f"Number of choices for conductor 9: {choices_9}")

    # For N = 7
    print("\n--- Component for conductor 7 ---")
    group_structure_7 = "C_6"
    num_prim_7 = phi(7) - 1
    print(f"The group of characters mod 7 is isomorphic to {group_structure_7}. It has {phi(7)} members.")
    print(f"A character is primitive mod 7 if it's not the principal character (conductor 1).")
    print(f"Number of primitive characters is phi(7) - 1 = {phi(7)} - 1 = {num_prim_7}.")
    print("The orders of characters mod 7 divide 6. All 5 primitive characters have orders (2, 3, 6) that divide 6.")
    choices_7 = num_prim_7
    print(f"Number of choices for conductor 7: {choices_7}")
    
    # For N = 11
    print("\n--- Component for conductor 11 ---")
    group_structure_11 = "C_10"
    print(f"The group of characters mod 11 is isomorphic to {group_structure_11}. It has {phi(11)} members.")
    print("A character must be primitive (order > 1) and its order must divide 6.")
    print("The order must divide gcd(10, 6) = 2. As it's primitive, the order must be 2.")
    choices_11 = phi(2)
    print(f"Number of characters of order 2 is phi(2) = {choices_11}.")
    print(f"Number of choices for conductor 11: {choices_11}")

    # For N = 13
    print("\n--- Component for conductor 13 ---")
    group_structure_13 = "C_12"
    print(f"The group of characters mod 13 is isomorphic to {group_structure_13}. It has {phi(13)} members.")
    print("A character must be primitive (order > 1) and its order must divide 6.")
    print("The order must divide gcd(12, 6) = 6. Possible orders are 2, 3, 6.")
    num_ord_2 = phi(2)
    num_ord_3 = phi(3)
    num_ord_6 = phi(6)
    choices_13 = num_ord_2 + num_ord_3 + num_ord_6
    print(f"Number of choices = (chars of order 2) + (chars of order 3) + (chars of order 6)")
    print(f"= phi(2) + phi(3) + phi(6) = {num_ord_2} + {num_ord_3} + {num_ord_6} = {choices_13}.")
    print(f"Number of choices for conductor 13: {choices_13}")

    print("\nStep 3: Verify overall order and calculate total")
    print("The order of the product character is the lcm of the orders of component characters:")
    print("lcm(ord(chi_4), ord(chi_9), ord(chi_7), ord(chi_11), ord(chi_13))")
    print("Possible orders are:")
    print("ord(chi_4) = 2")
    print("ord(chi_9) is 3 or 6")
    print("ord(chi_7) is 2, 3, or 6")
    print("ord(chi_11) = 2")
    print("ord(chi_13) is 2, 3, or 6")
    print("\nFor the lcm to be 6, it must be divisible by 2 and 3.")
    print("- It is always divisible by 2 because ord(chi_4)=2.")
    print("- It is always divisible by 3 because ord(chi_9) is always 3 or 6.")
    print("Since all component orders divide 6, the lcm is always exactly 6.")
    
    total_count = choices_4 * choices_9 * choices_7 * choices_11 * choices_13
    print("\nThus, any combination of the selected component characters results in a valid character of order 6.")
    print("The total number is the product of the number of choices for each component.")
    print(f"Total number = {choices_4} (for 4) * {choices_9} (for 9) * {choices_7} (for 7) * {choices_11} (for 11) * {choices_13} (for 13)")
    print(f"The final equation is: {choices_4} * {choices_9} * {choices_7} * {choices_11} * {choices_13} = {total_count}")

    return total_count

final_answer = solve_character_count()
print(f"<<<{final_answer}>>>")
