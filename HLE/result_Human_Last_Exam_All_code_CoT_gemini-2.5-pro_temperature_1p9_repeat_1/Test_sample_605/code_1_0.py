import math

# Step 1: Define weights and calculate the degree 'd'
weights = [22, 29, 49, 50, 75]
d = sum(weights)

# Step 2: Define the formula and calculate its components
# Formula: I = (2 * (d^3 - sum(w_i^3))) / (3 * d)
d_cubed = d**3
sum_w_cubed = sum([w**3 for w in weights])

# Calculate components for the final equation printout
w_cubed_list = [w**3 for w in weights]
numerator = 2 * (d_cubed - sum_w_cubed)
denominator = 3 * d

# Step 3: Output the equation with each number plugged in
print("The Crawley-Nordström invariant (I) is calculated as follows:")
print(f"I = (2 * ({d}^3 - ({' + '.join([str(w)+'^3' for w in weights])}))) / (3 * {d})")
print(f"I = (2 * ({d_cubed} - ({' + '.join([str(wc) for wc in w_cubed_list])}))) / ({denominator})")
print(f"I = (2 * ({d_cubed} - {sum_w_cubed})) / {denominator}")
print(f"I = (2 * {d_cubed - sum_w_cubed}) / {denominator}")
print(f"I = {numerator} / {denominator}")

# Step 4: Compute and print the final result
result = numerator / denominator
print(f"I = {result}")
