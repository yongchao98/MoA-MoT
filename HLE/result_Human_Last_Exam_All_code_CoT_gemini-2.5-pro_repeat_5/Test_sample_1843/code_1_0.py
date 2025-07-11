import math

def phi(n):
    """
    Computes Euler's totient function phi(n), which counts the number of
    positive integers up to a given integer n that are relatively prime to n.
    This is also the number of generators of a cyclic group of order n.
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
    Finds the number of primitive Dirichlet characters of conductor N and order k.
    """
    N = 36036
    order = 6

    print(f"We want to find the number of primitive Dirichlet characters of conductor N = {N} and order {order}.")
    print("Step 1: Factorize N into its prime power components.")
    print(f"N = {N} = 4 * 9 * 7 * 11 * 13\n")

    print("Step 2: A character chi modulo N is primitive if and only if it is a product of primitive characters")
    print("chi_m for each prime power factor m of N. The order of chi is the lcm of the orders of its components.\n")

    print(f"Step 3: We count the number of choices for each primitive component character chi_m such that its order divides {order}.")

    # For m = 4 = 2^2
    # The group of characters mod 4 is isomorphic to C_2. There is one primitive character of order 2.
    count_4 = 1 # phi(2) is not correct here, it's just 1.
    print(f"For m = 4: There is 1 primitive character. Its order is 2. Since 2 divides 6, we have {count_4} choice.")

    # For m = 9 = 3^2
    # Primitive characters mod 9 have orders that do not divide phi(3)=2. In C_6, these are orders 3 and 6.
    # Both 3 and 6 divide the target order 6.
    count_9 = phi(3) + phi(6)
    print(f"For m = 9: Primitive characters have orders 3 and 6. Both divide 6. Number of choices = phi(3) + phi(6) = {phi(3)} + {phi(6)} = {count_9}.")

    # For m = 7
    # Primitive characters mod 7 are all non-trivial characters. Orders are 2, 3, 6. All divide 6.
    count_7 = phi(2) + phi(3) + phi(6)
    print(f"For m = 7: Primitive characters can have orders 2, 3, 6. All divide 6. Number of choices = phi(2) + phi(3) + phi(6) = {phi(2)} + {phi(3)} + {phi(6)} = {count_7}.")

    # For m = 11
    # Primitive characters mod 11 have orders 2, 5, 10. Only order 2 divides 6.
    count_11 = phi(2)
    print(f"For m = 11: Primitive characters can have orders 2, 5, 10. Only order 2 divides 6. Number of choices = phi(2) = {count_11}.")

    # For m = 13
    # Primitive characters mod 13 have orders 2, 3, 4, 6, 12. Those dividing 6 are 2, 3, 6.
    count_13 = phi(2) + phi(3) + phi(6)
    print(f"For m = 13: Primitive characters can have orders 2, 3, 4, 6, 12. Those dividing 6 are 2, 3, 6. Number of choices = phi(2) + phi(3) + phi(6) = {phi(2)} + {phi(3)} + {phi(6)} = {count_13}.\n")

    print("Step 4: Verify that the lcm of the orders of any combination of these characters is exactly 6.")
    print("The order of chi_4 is 2, and the order of chi_11 is 2. So the lcm is always a multiple of 2.")
    print("The order of chi_9 is 3 or 6, so it's always a multiple of 3. The lcm is therefore a multiple of 3.")
    print("Thus, the final order is a multiple of lcm(2, 3) = 6.")
    print("Since all component orders divide 6, their lcm must also divide 6.")
    print("This means the order of any resulting character is exactly 6.\n")

    print("Step 5: The total number is the product of the number of choices for each component.")
    total = count_4 * count_9 * count_7 * count_11 * count_13
    print("The final calculation is:")
    print(f"{count_4} * {count_9} * {count_7} * {count_11} * {count_13} = {total}")

solve()