import math

def phi(n):
    """
    Calculates Euler's totient function phi(n).
    This function counts the number of positive integers up to a given integer n
    that are relatively prime to n.
    """
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

def count_characters(conductor, required_order_lcm):
    """
    Calculates the number of primitive Dirichlet characters for a given conductor and order.
    """
    N = conductor
    print(f"Finding the number of primitive Dirichlet characters of conductor N = {N} and order {required_order_lcm}.\n")

    # Step 1: Prime Factorization
    factors = {4: 2, 9: 2, 7: 1, 11: 1, 13: 1} # N = 36036 = 4 * 9 * 7 * 11 * 13
    print(f"Step 1: The prime factorization of N is 2^2 * 3^2 * 7 * 11 * 13.")
    print("This means we need to consider primitive characters modulo q = 4, 9, 7, 11, and 13.\n")

    # Step 2 & 3: Character and Order Decomposition
    print("Step 2: A character chi is primitive modulo N if and only if it is a product of primitive characters")
    print("chi = chi_4 * chi_9 * chi_7 * chi_11 * chi_13, where each chi_q is primitive modulo q.")
    print(f"The order of chi must be {required_order_lcm}, which means lcm(ord(chi_4), ..., ord(chi_13)) = {required_order_lcm}.\n")

    print("Step 3: Count the number of choices for each component character whose order divides 6.\n")

    # For q = 4 = 2^2
    # There is 1 primitive character modulo 4, and its order is 2.
    # Since 2 divides 6, this character is a valid choice.
    count_4 = 1
    print(f"For conductor q = 4: There is {count_4} primitive character. Its order is 2. (k_4 = 2)")

    # For q = 9 = 3^2
    # Primitive characters mod 9 have orders that are multiples of 3 and divide phi(9)=6.
    # Possible orders are 3 and 6. Both divide 6.
    # Number of characters of order 3 is phi(3)=2.
    # Number of characters of order 6 is phi(6)=2.
    count_9 = phi(3) + phi(6)
    print(f"For conductor q = 9: Primitive characters must have orders 3 or 6. There are {phi(3)} of order 3 and {phi(6)} of order 6. Total choices = {count_9}. (k_9 in {{3, 6}})")

    # For q = 7 (prime)
    # Primitive characters mod 7 have orders > 1 that divide phi(7)=6.
    # Possible orders are 2, 3, 6. All divide 6.
    # Number of characters of order 2 is phi(2)=1.
    # Number of characters of order 3 is phi(3)=2.
    # Number of characters of order 6 is phi(6)=2.
    count_7 = phi(2) + phi(3) + phi(6)
    print(f"For conductor q = 7: Primitive characters can have orders 2, 3, or 6. There are {phi(2)}+{phi(3)}+{phi(6)} = {count_7} choices. (k_7 in {{2, 3, 6}})")

    # For q = 11 (prime)
    # Primitive characters mod 11 have orders > 1 that divide phi(11)=10.
    # The order must also divide 6. The only common divisor > 1 is 2.
    # Number of characters of order 2 is phi(2)=1.
    count_11 = phi(2)
    print(f"For conductor q = 11: The primitive character's order must divide both 10 and 6. The only option > 1 is 2. There is {count_11} choice. (k_11 = 2)")

    # For q = 13 (prime)
    # Primitive characters mod 13 have orders > 1 that divide phi(13)=12.
    # The order must also divide 6. Possible orders are 2, 3, 6.
    # Number of characters of order 2 is phi(2)=1.
    # Number of characters of order 3 is phi(3)=2.
    # Number of characters of order 6 is phi(6)=2.
    count_13 = phi(2) + phi(3) + phi(6)
    print(f"For conductor q = 13: Primitive character orders can be 2, 3, or 6. There are {phi(2)}+{phi(3)}+{phi(6)} = {count_13} choices. (k_13 in {{2, 3, 6}})")
    print("-" * 20)

    # Step 4: Combine the counts
    print("Step 4: Combine the results.\n")
    print("The order of the combined character is the lcm of the orders of its components.")
    print("Let the orders be (k_4, k_9, k_7, k_11, k_13).")
    print("Since k_4=2 and k_9 is a multiple of 3, the lcm must be a multiple of lcm(2, 3) = 6.")
    print("Since all individual orders divide 6, their lcm must also divide 6.")
    print("Therefore, the lcm is always exactly 6 for any valid combination of component characters.\n")
    print("The total number of primitive characters of order 6 is the product of the number of choices for each component.")

    total_count = count_4 * count_9 * count_7 * count_11 * count_13

    print("\nFinal Calculation:")
    print(f"Number of characters = (choices for q=4) * (choices for q=9) * (choices for q=7) * (choices for q=11) * (choices for q=13)")
    print(f"Number of characters = {count_4} * {count_9} * {count_7} * {count_11} * {count_13} = {total_count}")
    return total_count

# --- Main Execution ---
if __name__ == "__main__":
    N = 36036
    ORDER = 6
    answer = count_characters(N, ORDER)
    print(f"\nThus, the number of primitive Dirichlet characters of conductor {N} and order {ORDER} is {answer}.")
