import math

def phi(n):
    """Computes Euler's totient function."""
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

# Step 4: Calculate the number of choices for each primitive component character
# whose order divides 6.

# For conductor 4: primitive characters have order 2. phi(2)=1.
choices_4 = phi(2)
print(f"For conductor 4, primitive characters must have order 2. Number of choices: phi(2) = {choices_4}")

# For conductor 9: primitive characters have order 3 or 6.
choices_9 = phi(3) + phi(6)
print(f"For conductor 9, primitive characters can have orders 3 or 6. Number of choices: phi(3) + phi(6) = {phi(3)} + {phi(6)} = {choices_9}")

# For conductor 7: primitive characters can have orders 2, 3, or 6.
choices_7 = phi(2) + phi(3) + phi(6)
print(f"For conductor 7, primitive characters can have orders 2, 3, or 6. Number of choices: phi(2) + phi(3) + phi(6) = {phi(2)} + {phi(3)} + {phi(6)} = {choices_7}")

# For conductor 11: primitive characters' orders must divide 6. Only order 2 is possible.
choices_11 = phi(2)
print(f"For conductor 11, primitive characters' order must divide 6. Only order 2 is possible. Number of choices: phi(2) = {choices_11}")

# For conductor 13: primitive characters' orders must divide 6. Possible orders are 2, 3, 6.
choices_13 = phi(2) + phi(3) + phi(6)
print(f"For conductor 13, primitive characters' orders must divide 6. Possible orders are 2, 3, 6. Number of choices: phi(2) + phi(3) + phi(6) = {phi(2)} + {phi(3)} + {phi(6)} = {choices_13}")

# Step 6: Final Calculation
total_chars = choices_4 * choices_9 * choices_7 * choices_11 * choices_13

print("\nThe total number of primitive Dirichlet characters of conductor 36036 and order 6 is the product of the number of choices for each component.")
print(f"Total number = {choices_4} * {choices_9} * {choices_7} * {choices_11} * {choices_13}")
print(f"Total number = {total_chars}")