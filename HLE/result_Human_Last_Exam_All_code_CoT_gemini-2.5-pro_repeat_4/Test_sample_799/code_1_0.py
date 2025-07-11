import math

# Step 1: Choose a large prime number p.
# This determines the size of our integer Sidon set.
p = 1009  # A reasonably large prime number

# Step 2: Define the properties of our constructed Sidon set of digits, D.
# Number theory allows us to construct an integer Sidon set D with p elements,
# where the maximum element is approximately p^2.
size_of_D = p
max_element_in_D = p**2 # An approximation for the largest element in D

# Step 3: Define the base 'q' for our number system.
# To ensure the resulting fractal is a Sidon set, we must prevent carries
# during the addition of two numbers. A base q > 2 * max(D) is sufficient.
# We will use q = 2 * p^2 for this calculation.
base_q = 2 * p**2

# Step 4: Calculate the Hausdorff dimension.
# The dimension of a set whose members are formed by digits from D in base q
# is given by the formula: log(|D|) / log(q).
hausdorff_dimension = math.log(size_of_D) / math.log(base_q)

# Step 5: Print the explanation and the result.
# As p gets larger, this value gets arbitrarily close to 0.5.
print("Illustrative Calculation for the Hausdorff Dimension of a Sidon Set")
print("---------------------------------------------------------------------")
print(f"We construct the set using a prime p = {p}.")
print(f"Size of the digit set |D| = p = {size_of_D}")
print(f"Base for the number system q = 2 * p^2 = {base_q}")
print("\nThe dimension is calculated using the formula: log(|D|) / log(q)")
print(f"\nlog({size_of_D}) / log({base_q}) = {hausdorff_dimension:.8f}")
print("\nThis demonstrates that we can construct Sidon sets with a dimension very close to 0.5.")
print("The theoretical maximum dimension is exactly 0.5.")
