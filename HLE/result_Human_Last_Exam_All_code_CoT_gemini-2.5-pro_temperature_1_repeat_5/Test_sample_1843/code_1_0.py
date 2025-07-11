import math

def phi(n):
    """
    Calculates Euler's totient function phi(n).
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer")
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
    Calculates the number of primitive Dirichlet characters of a given
    conductor and order.
    """
    N = 36036
    ORDER_WANTED = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor N = {N} and order {ORDER_WANTED}.")
    print("-" * 50)

    # Step 1: Factor N
    # N = 36036 = 36 * 1001 = (2^2 * 3^2) * (7 * 11 * 13)
    prime_powers = [4, 9, 7, 11, 13]
    
    print(f"The conductor N factorizes into prime powers: {prime_powers}")
    print("A character is primitive of conductor N if and only if it is a product of primitive characters with conductors from this set.")
    print("The order of the character is the lcm of the orders of its components.")
    print("-" * 50)

    # Step 2: For each component, find the number of primitive characters whose order divides 6.
    print(f"We analyze each component to find primitive characters whose order divides {ORDER_WANTED}:")
    
    counts = {}

    # Component q=4 (2^2)
    # The group of characters mod 4 is cyclic of order phi(4)=2.
    # The number of primitive characters is phi(4) - phi(2) = 1. Its order is 2.
    # Since 2 divides 6, this character is a valid component.
    c4 = 1
    counts[4] = c4
    print(f"\nFor conductor q=4:")
    print(f"  There is {c4} primitive character. Its order is 2, which divides 6.")

    # Component q=9 (3^2)
    # The group of characters mod 9 is cyclic of order phi(9)=6.
    # The primitive characters are those not induced from mod 3, i.e., those with orders
    # that are not 1 or 2. The possible orders are 3 and 6. Both divide 6.
    c9_ord3 = phi(3)
    c9_ord6 = phi(6)
    c9 = c9_ord3 + c9_ord6
    counts[9] = c9
    print(f"\nFor conductor q=9:")
    print("  Primitive characters have orders 3 or 6. Both divide 6.")
    print(f"  - Number of characters of order 3: {c9_ord3}")
    print(f"  - Number of characters of order 6: {c9_ord6}")
    print(f"  Total valid choices for q=9: {c9}")

    # Component q=7
    # The group of characters mod 7 is cyclic of order phi(7)=6.
    # All non-principal characters are primitive. Their orders are 2, 3, 6. All divide 6.
    c7_ord2 = phi(2)
    c7_ord3 = phi(3)
    c7_ord6 = phi(6)
    c7 = c7_ord2 + c7_ord3 + c7_ord6
    counts[7] = c7
    print(f"\nFor conductor q=7:")
    print("  Primitive characters have orders 2, 3, or 6. All divide 6.")
    print(f"  - Number of characters of order 2: {c7_ord2}")
    print(f"  - Number of characters of order 3: {c7_ord3}")
    print(f"  - Number of characters of order 6: {c7_ord6}")
    print(f"  Total valid choices for q=7: {c7}")

    # Component q=11
    # The group of characters mod 11 is cyclic of order phi(11)=10.
    # Primitive character orders are 2, 5, 10. Only order 2 divides 6.
    c11_ord2 = phi(2)
    c11 = c11_ord2
    counts[11] = c11
    print(f"\nFor conductor q=11:")
    print("  Primitive character orders are 2, 5, 10. Only order 2 divides 6.")
    print(f"  - Number of characters of order 2: {c11}")
    print(f"  Total valid choices for q=11: {c11}")

    # Component q=13
    # The group of characters mod 13 is cyclic of order phi(13)=12.
    # Primitive character orders are 2, 3, 4, 6, 12. Orders 2, 3, 6 divide 6.
    c13_ord2 = phi(2)
    c13_ord3 = phi(3)
    c13_ord6 = phi(6)
    c13 = c13_ord2 + c13_ord3 + c13_ord6
    counts[13] = c13
    print(f"\nFor conductor q=13:")
    print("  Primitive character orders are 2, 3, 4, 6, 12. Orders 2, 3, and 6 divide 6.")
    print(f"  - Number of characters of order 2: {c13_ord2}")
    print(f"  - Number of characters of order 3: {c13_ord3}")
    print(f"  - Number of characters of order 6: {c13_ord6}")
    print(f"  Total valid choices for q=13: {c13}")

    print("-" * 50)
    # Step 3 & 4: Combine counts and apply lcm condition.
    total_combinations = 1
    for q in prime_powers:
        total_combinations *= counts[q]

    print("The total number of combinations of component characters whose orders all divide 6 is the product of the individual counts:")
    equation = " * ".join(str(c) for c in counts.values())
    print(f"Total combinations = {equation} = {total_combinations}")
    print()

    print("For the final character to have order 6, the lcm of the component orders must be exactly 6.")
    print("This means the lcm must be divisible by 2 and 3.")
    print("\n- Divisibility by 2:")
    print("  The component for conductor 4 must have order 2. Thus, the lcm is always divisible by 2.")
    print("\n- Divisibility by 3:")
    print("  The component for conductor 9 must have order 3 or 6. Thus, the lcm is always divisible by 3.")
    print()
    
    print("Since any combination of our selected components will result in an lcm that is a multiple of 2 and 3, and since all component orders divide 6, the final lcm must be exactly 6.")
    print(f"Therefore, the final answer is the total number of combinations.")
    
    final_answer = total_combinations
    print("-" * 50)
    print(f"The number of primitive Dirichlet characters of conductor {N} and order {ORDER_WANTED} is {final_answer}.")
    return final_answer

if __name__ == '__main__':
    solve()