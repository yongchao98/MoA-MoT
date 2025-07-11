import math

def phi(n):
    """
    Calculates Euler's totient function.
    Not strictly required as we only need phi(2), phi(3), phi(6),
    but included for completeness.
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

def solve():
    """
    Finds the number of primitive Dirichlet characters of conductor N=36036 and order 6.
    """
    N = 36036
    # Step 1: Prime factorization of N
    # N = 36036 = 36 * 1001 = 4 * 9 * 7 * 11 * 13 = 2^2 * 3^2 * 7 * 11 * 13
    print(f"The conductor is N = 36036.")
    print("Its prime factorization is N = 2^2 * 3^2 * 7 * 11 * 13.\n")

    print("A primitive character chi mod N can be decomposed into primitive characters for each prime power factor:")
    print("chi = chi_4 * chi_9 * chi_7 * chi_11 * chi_13")
    print("The order of chi must be lcm(ord(chi_4), ord(chi_9), ...) = 6.\n")

    print("We count the number of choices for each component character whose order divides 6.")

    # Step 2: Count choices for each prime power factor.

    # Modulus 4 (2^2)
    # The character group mod 4, (Z/4Z)*, is cyclic of order phi(4) = 2.
    # It has one primitive character of order 2. Its order (2) divides 6.
    c_4 = 1
    print(f"For modulus 4 (2^2): There is {c_4} primitive character of order 2.")

    # Modulus 9 (3^2)
    # The character group mod 9, (Z/9Z)*, is cyclic of order phi(9) = 6.
    # Primitive characters mod 9 have orders 3 and 6. Both divide 6.
    # Number of primitive characters of order 3 = phi(3) = 2.
    # Number of primitive characters of order 6 = phi(6) = 2.
    c_9 = phi(3) + phi(6)
    print(f"For modulus 9 (3^2): There are {phi(3)} primitive characters of order 3 and {phi(6)} of order 6. Total = {c_9} choices.")

    # Modulus 7
    # The character group mod 7 is cyclic of order phi(7) = 6.
    # All non-principal characters are primitive. We need those with order dividing 6.
    # Orders are 2, 3, 6.
    # Number of primitive characters of order 2 = phi(2) = 1.
    # Number of primitive characters of order 3 = phi(3) = 2.
    # Number of primitive characters of order 6 = phi(6) = 2.
    c_7 = phi(2) + phi(3) + phi(6)
    print(f"For modulus 7: There are {phi(2)} primitive character of order 2, {phi(3)} of order 3, and {phi(6)} of order 6. Total = {c_7} choices.")

    # Modulus 11
    # The character group mod 11 is cyclic of order phi(11) = 10.
    # We need primitive characters whose order divides 6.
    # The divisors of 10 are {1, 2, 5, 10}. The only relevant order > 1 is 2.
    # Number of primitive characters of order 2 = phi(2) = 1.
    c_11 = phi(2)
    print(f"For modulus 11: There is {c_11} primitive character of order 2.")

    # Modulus 13
    # The character group mod 13 is cyclic of order phi(12) = 12.
    # We need primitive characters whose order divides 6. Relevant orders are 2, 3, 6.
    # Number of primitive characters of order 2 = phi(2) = 1.
    # Number of primitive characters of order 3 = phi(3) = 2.
    # Number of primitive characters of order 6 = phi(6) = 2.
    c_13 = phi(2) + phi(3) + phi(6)
    print(f"For modulus 13: There are {phi(2)} primitive character of order 2, {phi(3)} of order 3, and {phi(6)} of order 6. Total = {c_13} choices.\n")

    # Step 3: Combine the results.
    print("To get a character of order 6, the lcm of the component orders must be 6.")
    print("- The character for mod 4 must have order 2. So the lcm is always divisible by 2.")
    print("- The character for mod 9 must have order 3 or 6. So the lcm is always divisible by 3.")
    print("- Since all component orders divide 6, their lcm must divide 6.")
    print("Any number that is a multiple of 2, a multiple of 3, and divides 6 must be 6.")
    print("Thus, any combination of the selected component characters will have a resulting order of 6.\n")

    total_count = c_4 * c_9 * c_7 * c_11 * c_13
    
    print("The total number of primitive Dirichlet characters of conductor 36036 and order 6 is the product of the number of choices for each factor:")
    print(f"{c_4} * {c_9} * {c_7} * {c_11} * {c_13} = {total_count}")

solve()