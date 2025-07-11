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

# Step 1: Count the number of choices for each primitive component character
# such that its order divides 6.

# For conductor 4 (2^2): 1 primitive character of order 2. 2|6.
num_chi_4 = 1

# For conductor 9 (3^2): primitive characters have orders 3 or 6. Both divide 6.
# Number of choices = phi(3) + phi(6)
num_chi_9 = phi(3) + phi(6)

# For conductor 7: primitive characters have orders 2, 3, or 6. All divide 6.
# Number of choices = phi(2) + phi(3) + phi(6)
num_chi_7 = phi(2) + phi(3) + phi(6)

# For conductor 11: primitive characters have orders 2, 5, or 10. Only 2 divides 6.
# Number of choices = phi(2)
num_chi_11 = phi(2)

# For conductor 13: primitive characters have orders 2, 3, 4, 6, 12.
# Those dividing 6 are 2, 3, and 6.
# Number of choices = phi(2) + phi(3) + phi(6)
num_chi_13 = phi(2) + phi(3) + phi(6)

# Step 2: Multiply the number of choices for each component.
# As explained in the text, the lcm is guaranteed to be 6.
total_characters = num_chi_4 * num_chi_9 * num_chi_7 * num_chi_11 * num_chi_13

print("The number of primitive Dirichlet characters of conductor 36036 and order 6 is the product of the number of choices for each primitive component:")
print(f"Choices for N=4: {num_chi_4}")
print(f"Choices for N=9: {phi(3)} + {phi(6)} = {num_chi_9}")
print(f"Choices for N=7: {phi(2)} + {phi(3)} + {phi(6)} = {num_chi_7}")
print(f"Choices for N=11: {phi(2)} = {num_chi_11}")
print(f"Choices for N=13: {phi(2)} + {phi(3)} + {phi(6)} = {num_chi_13}")
print("\nFinal Calculation:")
print(f"{num_chi_4} * {num_chi_9} * {num_chi_7} * {num_chi_11} * {num_chi_13} = {total_characters}")
