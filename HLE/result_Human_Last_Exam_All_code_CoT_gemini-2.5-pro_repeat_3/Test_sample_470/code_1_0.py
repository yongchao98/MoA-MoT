import math

# --- Problem Parameters ---
# Let D = (C_2)^5 be the defect group.
# This is an elementary abelian 2-group, which can be viewed as a vector space
# of dimension 5 over the field with 2 elements, F_2.
dim_D = 5

# Let E be the inertial quotient, with |E|=5.
order_E = 5

# The field F has characteristic 2.
p = 2

print("Step 1: Calculate l(B), the number of irreducible Brauer characters.")
# The number of irreducible Brauer characters, l(B), is the number of irreducible
# modules for the inertial quotient algebra F[E].
# The characteristic of F (p=2) does not divide the order of E (|E|=5).
# Thus, the algebra F[E] is semisimple.
# For a semisimple algebra over a splitting field F, the number of simple modules
# equals the number of conjugacy classes of the group E.
# Since E has prime order 5, it is cyclic and therefore abelian. The number of
# conjugacy classes is simply its order.
l_B = order_E
print(f"The group E is abelian, so the number of its conjugacy classes is |E|.")
print(f"Therefore, l(B) = |E| = {l_B}.\n")


print("Step 2: Calculate k(B), the number of irreducible ordinary characters.")
# For a block with an abelian defect group D, k(B) is the number of orbits of
# the inertial quotient E acting on the set of irreducible characters of D, Irr(D).
# We use Burnside's Lemma to count these orbits:
# Number of orbits = (1/|E|) * sum_{g in E} |Irr(D)^g|
# where |Irr(D)^g| is the number of characters in Irr(D) fixed by the element g.

# The character group Irr(D) is isomorphic to D, so |Irr(D)| = |D| = 2^dim_D.
size_Irr_D = 2**dim_D

# We need to find the number of fixed points for each element of E.
# The action of E on Irr(D) is dual to its action on D. The number of fixed
# points for an element g is the same for both actions.

# Case 1: The identity element 'e' of E.
# The identity element fixes all characters.
fixed_points_identity = size_Irr_D
num_identity = 1
print(f"The identity element of E fixes all {fixed_points_identity} characters of D.")

# Case 2: A non-identity element 'g' of E.
# Since |E|=5, any non-identity element g has order 5. Its action on D corresponds
# to an element of order 5 in Aut(D) = GL(5, 2).
# The number of points in D fixed by g is the size of the kernel of (A - I),
# where A is the matrix for g. This is the eigenspace for eigenvalue 1.
# The minimal polynomial of A divides x^5-1 = (x+1)(x^4+x^3+x^2+x+1) over F_2.
# Since A is not the identity, its minimal polynomial is the irreducible x^4+x^3+x^2+x+1.
# The characteristic polynomial of A must be (x+1)(x^4+x^3+x^2+x+1).
# The dimension of the eigenspace for eigenvalue 1 is the multiplicity of the root 1, which is 1.
# So, the subspace of fixed points has dimension 1, containing 2^1 = 2 elements.
fixed_points_non_identity = 2**1
num_non_identity = order_E - 1
print(f"Any non-identity element of E fixes {fixed_points_non_identity} characters of D.")

# Apply Burnside's Lemma
sum_of_fixed_points = (num_identity * fixed_points_identity) + (num_non_identity * fixed_points_non_identity)
k_B = sum_of_fixed_points // order_E

print(f"\nUsing Burnside's Lemma, k(B) is the average number of fixed points:")
print(f"k(B) = (1/{order_E}) * ({num_identity} * {fixed_points_identity} + {num_non_identity} * {fixed_points_non_identity})")
print(f"k(B) = (1/{order_E}) * ({sum_of_fixed_points}) = {k_B}\n")


print("Step 3: Compute the final result k(B) - l(B).")
result = k_B - l_B
print(f"The value of k(B) - l(B) is: {k_B} - {l_B} = {result}")
