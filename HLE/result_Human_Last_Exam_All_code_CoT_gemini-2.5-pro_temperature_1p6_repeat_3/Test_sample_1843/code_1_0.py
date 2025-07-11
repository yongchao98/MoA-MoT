import math

def euler_phi(n):
    """Calculates Euler's totient function."""
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

# The number of primitive characters for each prime power factor p^k whose order divides 6.
# Calculated as Sum(phi(d)) for valid orders d.

# For conductor 4 (2^2): primitive order is 2. phi(2) = 1
count_4 = euler_phi(2)

# For conductor 9 (3^2): primitive orders are 3, 6. phi(3) + phi(6) = 2 + 2 = 4
count_9 = euler_phi(3) + euler_phi(6)

# For conductor 7: primitive orders are 2, 3, 6. phi(2) + phi(3) + phi(6) = 1 + 2 + 2 = 5
count_7 = euler_phi(2) + euler_phi(3) + euler_phi(6)

# For conductor 11: primitive order must divide gcd(10,6)=2. order is 2. phi(2) = 1
count_11 = euler_phi(2)

# For conductor 13: primitive order must divide gcd(12,6)=6. orders are 2, 3, 6. phi(2)+phi(3)+phi(6) = 1+2+2=5
count_13 = euler_phi(2) + euler_phi(3) + euler_phi(6)

component_counts = [count_4, count_9, count_7, count_11, count_13]

# Calculate the total number by multiplying the counts for each component.
total_count = math.prod(component_counts)

# Format the final equation string
equation_str = " * ".join(map(str, component_counts))

print("The number of primitive Dirichlet characters of conductor N=36036 and order 6 is the product of the number of choices for each prime power component.")
print(f"The number of choices for conductors 4, 9, 7, 11, and 13 are {component_counts[0]}, {component_counts[1]}, {component_counts[2]}, {component_counts[3]}, and {component_counts[4]} respectively.")
print("\nThe final calculation is:")
print(f"{equation_str} = {total_count}")
