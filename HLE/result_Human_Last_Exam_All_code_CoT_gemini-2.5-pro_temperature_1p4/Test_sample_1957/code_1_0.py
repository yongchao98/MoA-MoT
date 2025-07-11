import math

# Step 1 & 2: Determine the degree D of the polynomial P.
# Based on the geometric analysis, the condition is equivalent to an algebraic
# equation P=0 where the polynomial P has a total degree of D=8.
D = 8

# Step 3: Minimize 1000m + n subject to m * n = D.
# We are looking for positive integers m and n.
# We can list all possible integer pairs (m, n) such that m * n = 8.
factors = []
for i in range(1, int(math.sqrt(D)) + 1):
    if D % i == 0:
        factors.append((i, D // i))
        if i*i != D:
            factors.append((D // i, i))

# Sort the factors by m for clarity.
factors.sort()

min_val = float('inf')
best_pair = (None, None)

print(f"The degree of the defining polynomial is D = {D}.")
print("We need to minimize 1000m + n subject to m * n = D.")
print("The possible integer pairs (m, n) are:")

for m, n in factors:
    value = 1000 * m + n
    print(f"  m = {m}, n = {n}: 1000 * {m} + {n} = {value}")
    if value < min_val:
        min_val = value
        best_pair = (m, n)

print(f"\nThe minimum value is obtained for the pair (m, n) = {best_pair}.")
print(f"The minimum value of 1000m + n is {min_val}.")