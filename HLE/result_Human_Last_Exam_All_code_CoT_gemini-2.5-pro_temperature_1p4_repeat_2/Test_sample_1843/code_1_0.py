import math

def phi(n):
    """Computes Euler's totient function, which counts the number of elements of a given order in a cyclic group."""
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

print("Let N be the number of primitive Dirichlet characters of conductor 36036 and order 6.")
print("The conductor is N = 36036 = 2^2 * 3^2 * 7 * 11 * 13.")
print("A character chi is primitive mod N if it is a product of primitive characters for each prime power factor.")
print("The order of chi is the lcm of the orders of its components.\n")

# For each modulus, count primitive characters with orders dividing 6.
# c_m[k] is the number of primitive characters mod m of order k.

# Mod 4: Primitive characters have order 2.
c4 = {2: 1}
# Mod 9: Primitive characters have orders 3 and 6.
c9 = {3: phi(3), 6: phi(6)} # {3: 2, 6: 2}
# Mod 7: Primitive characters have orders 2, 3, 6.
c7 = {2: phi(2), 3: phi(3), 6: phi(6)} # {2: 1, 3: 2, 6: 2}
# Mod 11: Primitive characters have orders 2, 5, 10.
c11 = {2: phi(2), 5: phi(5), 10: phi(10)} # {2: 1, 5: 4, 10: 4}
# Mod 13: Primitive characters have orders 2, 3, 4, 6, 12.
c13 = {2: phi(2), 3: phi(3), 4: phi(4), 6: phi(6), 12: phi(12)} # {2: 1, 3: 2, 4: 2, 6: 2, 12: 4}

print("We count the number of choices for each component character (chi_4, chi_9, chi_7, chi_11, chi_13):")

# Number of combinations where all orders divide 6 (lcm divides 6)
n4_div_6 = c4.get(2, 0)
n9_div_6 = c9.get(3, 0) + c9.get(6, 0)
n7_div_6 = c7.get(2, 0) + c7.get(3, 0) + c7.get(6, 0)
n11_div_6 = c11.get(2, 0)
n13_div_6 = c13.get(2, 0) + c13.get(3, 0) + c13.get(6, 0)
C6 = n4_div_6 * n9_div_6 * n7_div_6 * n11_div_6 * n13_div_6

print(f"Number of choices where order divides 6 for mod 4: {n4_div_6}")
print(f"Number of choices where order divides 6 for mod 9: {n9_div_6}")
print(f"Number of choices where order divides 6 for mod 7: {n7_div_6}")
print(f"Number of choices where order divides 6 for mod 11: {n11_div_6}")
print(f"Number of choices where order divides 6 for mod 13: {n13_div_6}")
print(f"Total combinations where lcm of orders divides 6 (C_6) = {n4_div_6} * {n9_div_6} * {n7_div_6} * {n11_div_6} * {n13_div_6} = {C6}\n")

# Number of combinations where all orders divide 2 (lcm divides 2)
n4_div_2 = c4.get(2, 0)
n9_div_2 = c9.get(2, 0)
n7_div_2 = c7.get(2, 0)
n11_div_2 = c11.get(2, 0)
n13_div_2 = c13.get(2, 0)
C2 = n4_div_2 * n9_div_2 * n7_div_2 * n11_div_2 * n13_div_2
print(f"Total combinations where lcm of orders divides 2 (C_2) = {n4_div_2} * {n9_div_2} * {n7_div_2} * {n11_div_2} * {n13_div_2} = {C2}")

# Number of combinations where all orders divide 3 (lcm divides 3)
n4_div_3 = c4.get(3, 0)
n9_div_3 = c9.get(3, 0)
n7_div_3 = c7.get(3, 0)
n11_div_3 = c11.get(3, 0)
n13_div_3 = c13.get(3, 0)
C3 = n4_div_3 * n9_div_3 * n7_div_3 * n11_div_3 * n13_div_3
print(f"Total combinations where lcm of orders divides 3 (C_3) = {n4_div_3} * {n9_div_3} * {n7_div_3} * {n11_div_3} * {n13_div_3} = {C3}")

# Number of combinations where all orders divide 1 (lcm divides 1)
C1 = 0
print(f"Total combinations where lcm of orders divides 1 (C_1) = {C1}\n")

# By the Principle of Inclusion-Exclusion, the number of characters with order exactly 6 is:
# N = C_6 - C_2 - C_3 + C_1
result = C6 - C2 - C3 + C1

print("By the Principle of Inclusion-Exclusion, the final number is:")
print(f"N = C_6 - C_2 - C_3 + C_1 = {C6} - {C2} - {C3} + {C1} = {result}")
