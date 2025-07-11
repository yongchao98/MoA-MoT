# Step 1: Define the given values as variables.
# l1 is the value assumed to be x11.
# c1 is the given constant.
l1 = (1 + 10**5)**5
c1 = 10**4

# Step 2: Calculate the numerator (x11 - 1).
numerator = l1 - 1

# Step 3: Calculate u1 by performing the integer division.
# Python's integers have arbitrary precision, avoiding overflow issues.
u1 = numerator // c1

# Step 4: Print the final equation with all the numbers.
print(f"From the matrix equation, we derive the formula for u1: u1 = (x11 - 1) / c1")
print(f"Assuming x11 = l1, and with the given values:")
print(f"l1 = (1 + 10^5)^5 = {l1}")
print(f"c1 = 10^4 = {c1}")
print(f"The equation for u1 is:")
print(f"({l1} - 1) / {c1} = {u1}")