import math

def euler_totient(n):
    """
    Calculates Euler's totient function phi(n).
    It counts the numbers from 1 to n that are relatively prime to n.
    """
    count = 0
    for i in range(1, n + 1):
        if math.gcd(i, n) == 1:
            count += 1
    return count

# Define the parameters for the reflection group of type H3
group_order = 120
coxeter_number = 10

# The problem is equivalent to finding the number of regular elements in H3.
# The formula is: N = |W| * phi(h) / h

# Calculate Euler's totient function for the Coxeter number
phi_h = euler_totient(coxeter_number)

# Calculate the total number of elements
num_elements = (group_order * phi_h) // coxeter_number

# Print the final equation and the result
print("The number of elements is found using the formula: |W| * phi(h) / h")
print(f"For the H3 group, the order |W| = {group_order}.")
print(f"The Coxeter number h = {coxeter_number}.")
print(f"Euler's totient function, phi({coxeter_number}), is {phi_h}.")
print(f"The final calculation is: {group_order} * {phi_h} / {coxeter_number} = {num_elements}")