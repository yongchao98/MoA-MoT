import math

# Step 1: Define the given parameters
p = 2
D_order = 2**5
E_order = 5

# Step 2: Calculate k(B), the number of ordinary characters.
# This is the number of E-orbits on Irr(D).
# We use Burnside's Lemma: (1/|E|) * sum_{g in E} |Fix(g)|
# |Fix(g)| = |C_D(g)|, the number of elements in D fixed by g.
# For g = identity (1 element): |C_D(g)| = |D| = 32
# For g != identity (4 elements): g fixes the 1D trivial submodule, so |C_D(g)| = 2^1 = 2
num_identity = 1
num_non_identity = E_order - 1

fixed_by_id = D_order
fixed_by_non_id = 2 # From the 1-dimensional trivial submodule

sum_of_fixed_points = (num_identity * fixed_by_id) + (num_non_identity * fixed_by_non_id)
k_B = sum_of_fixed_points // E_order

print(f"The number of ordinary characters k(B) is calculated as:")
print(f"k(B) = (1 * {fixed_by_id} + {num_non_identity} * {fixed_by_non_id}) / {E_order}")
print(f"k(B) = ({sum_of_fixed_points}) / {E_order}")
print(f"k(B) = {k_B}\n")

# Step 3: Calculate l(B), the number of Brauer characters.
# This is the number of irreducible F-characters of E.
# Since char(F) = 2 and |E| = 5, E is a 2'-group.
# Since F is large enough, it is a splitting field for E.
# For an abelian group over a splitting field, the number of irreducible characters is |E|.
l_B = E_order
print(f"The number of Brauer characters l(B) is the order of the inertial quotient E.")
print(f"l(B) = {l_B}\n")

# Step 4: Compute the difference k(B) - l(B)
difference = k_B - l_B

print("The value of k(B) - l(B) is:")
print(f"{k_B} - {l_B} = {difference}")
