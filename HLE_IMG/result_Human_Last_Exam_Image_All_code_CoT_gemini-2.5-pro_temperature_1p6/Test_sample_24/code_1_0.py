import numpy as np

# Step 1: Identify the manifold from its Heegaard diagram
# The diagram is a known representation of the Poincaré homology sphere.
manifold_name = "Poincaré homology sphere"
print(f"The three-manifold represented by the Heegaard diagram is the {manifold_name}.")
print("-" * 40)

# Step 2: Provide the fundamental group presentation
# The fundamental group of this manifold is the binary icosahedral group.
group_presentation = "<x, y | x³ = y⁵ = (xy)²>"
print(f"Its fundamental group, π₁(M), has the presentation:")
print(group_presentation)
print("-" * 40)

# Step 3: Calculate the first homology group H₁(M, ℤ)
# H₁(M, ℤ) is the abelianization of π₁(M). We derive its presentation matrix
# from the group relators.
# In an abelian group, the relations become (in additive notation):
# 1) x³ = y⁵   => 3x - 5y = 0
# 2) y⁵ = (xy)² = x²y² => y³ = x² => 2x - 3y = 0
# The coefficients form the presentation matrix.
print("To confirm it's a homology sphere, we calculate its first homology group, H₁(M, ℤ).")
print("The presentation matrix A for H₁(M, ℤ) is:")
A = np.array([[3, -5],
              [2, -3]])
print(A)
print("\nThe order of the group H₁(M, ℤ) is the absolute value of the determinant of A.")
print("-" * 40)

# Step 4: Show the determinant calculation
# The final equation demonstrates that the determinant is 1.
a, b = A[0, 0], A[0, 1]
c, d = A[1, 0], A[1, 1]
determinant = (a * d) - (b * c)

print("Calculation of the determinant:")
print(f"det(A) = ({a} * {d}) - ({b} * {c})")
print(f"       = {a * d} - ({b * c})")
print(f"       = {a * d} + {-b * c}")
print(f"       = {int(determinant)}")
print("-" * 40)

# Step 5: Conclude based on the result
if abs(int(determinant)) == 1:
    print("Since the determinant is 1, H₁(M, ℤ) is the trivial group.")
    print("This confirms the manifold is a homology sphere.")
else:
    print("The calculation shows the homology group is not trivial.")
