import math

def phi(n):
    """
    Computes Euler's totient function.
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
    Solves the problem of finding the number of primitive Dirichlet characters
    of conductor N=36036 and order 6.
    """
    N = 36036
    print(f"Finding the number of primitive Dirichlet characters of conductor N = {N} and order 6.")
    
    print("\nStep 1: Prime factorization of N.")
    print(f"N = 36036 = 2^2 * 3^2 * 7 * 11 * 13")
    
    print("\nStep 2: A primitive character chi mod N is a product of primitive characters mod p^a.")
    print("chi = chi_4 * chi_9 * chi_7 * chi_11 * chi_13")
    
    print("\nStep 3: The order of chi is the least common multiple (lcm) of the orders of its components.")
    print("We need lcm(ord(chi_4), ord(chi_9), ...) = 6.")
    print("This means each component's order must divide 6.")
    
    print("\nStep 4: Count the number of choices for each primitive component whose order divides 6.")
    
    # Modulo 4
    # The group of characters is isomorphic to (Z/4Z)*, which is cyclic of order 2.
    # It has one primitive character of order 2. 2 divides 6.
    count_4 = 1
    print(f"For modulus 4 (2^2): There is 1 primitive character (order 2). Its order divides 6. Count = {count_4}")
    
    # Modulo 9
    # The number of primitive characters is phi(9) - phi(3) = 6 - 2 = 4.
    # Their orders are 3 and 6, both divide 6.
    count_9 = phi(9) - phi(3)
    print(f"For modulus 9 (3^2): There are {phi(9)}-{phi(3)} = {count_9} primitive characters (orders 3, 6). All their orders divide 6. Count = {count_9}")
    
    # Modulo 7
    # Primitive characters are non-trivial. Orders are divisors of phi(7)=6, i.e. 2,3,6. All divide 6.
    count_7 = phi(7) - 1
    print(f"For modulus 7: There are {phi(7)}-1 = {count_7} primitive characters (orders 2, 3, 6). All their orders divide 6. Count = {count_7}")
    
    # Modulo 11
    # Primitive characters have orders dividing phi(11)=10. We need those that divide 6. Only order 2 works.
    # Number of characters of order 2 is phi(2)=1.
    count_11 = phi(2)
    print(f"For modulus 11: Primitive character orders must divide 10. Only order 2 also divides 6. There is phi(2) = {count_11} such character. Count = {count_11}")
    
    # Modulo 13
    # Primitive characters have orders dividing phi(13)=12. We need those that divide 6: orders 2, 3, 6.
    # Number of such characters is phi(2) + phi(3) + phi(6)
    count_13 = phi(2) + phi(3) + phi(6)
    print(f"For modulus 13: Primitive character orders must divide 12. Those dividing 6 are 2, 3, 6. Number of such characters is phi(2)+phi(3)+phi(6) = {phi(2)}+{phi(2)}+{phi(6)} = {count_13}. Count = {count_13}")

    print("\nStep 5: The total number of primitive characters with order dividing 6 is the product of these counts.")
    total_count = count_4 * count_9 * count_7 * count_11 * count_13
    print(f"Total count = {count_4} * {count_9} * {count_7} * {count_11} * {count_13} = {total_count}")
    
    print("\nStep 6: Final verification.")
    print("The order of any character chi_4 must be 2, so lcm is a multiple of 2.")
    print("The order of any character chi_9 must be 3 or 6, so lcm is a multiple of 3.")
    print("Any character counted has an order that is a multiple of lcm(2,3)=6 and also a divisor of 6.")
    print("Therefore, the order must be exactly 6.")
    
    print("\nThe number of primitive Dirichlet characters of conductor 36036 and order 6 is the total count.")
    print(f"\nFinal Answer: {total_count}")

solve()