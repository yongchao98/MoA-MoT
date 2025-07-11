# The fundamental group of the Poincaré sphere is SL(2, 5).
# We calculate its order using the formula for SL(2, F_q) where q=5.
q = 5

# The formula for the order is q * (q^2 - 1).
order = q * (q**2 - 1)

# Print the equation with the numbers substituted
print(f"The order of the fundamental group is calculated as:")
print(f"{q} * ({q}^2 - 1) = {q} * ({q*q} - 1) = {q} * ({q*q - 1}) = {order}")
print(f"\nThe fundamental group has {order} elements.")

# Another way to calculate this is to iterate through all matrices and count.
# This serves as a verification of the formula.
count = 0
for a in range(q):
    for b in range(q):
        for c in range(q):
            for d in range(q):
                # Check if the determinant ad - bc is 1 (mod q)
                if (a * d - b * c) % q == 1:
                    count += 1
print(f"Verification by direct counting: {count}")