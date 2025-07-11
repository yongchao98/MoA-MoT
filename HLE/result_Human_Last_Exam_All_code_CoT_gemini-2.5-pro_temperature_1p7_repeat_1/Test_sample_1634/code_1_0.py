# The problem asks for the smallest non-negative integer n such that
# there exists an n-point topological space that is not irreducible.

# A topological space X is NOT irreducible if it can be written as a
# finite union of proper closed subsets.
# X = Z1 U Z2 U ... U Zk, where each Zi is a closed subset and Zi != X.

# For n=0, the space is X = {}. It has no proper subsets, so it's irreducible.
# For n=1, the space is X = {p}. The only proper closed subset is {}.
# Any union of {} is still {}, which is not X. So, any 1-point space is irreducible.

# Let's consider n=2.
# Let the space X be a set with 2 points.
X = {1, 2}
n = len(X)
print(f"Let's test if a space with n = {n} points can be not irreducible.")
print(f"Consider the space X = {X}")

# We need to define a topology on X. Let's use the discrete topology,
# where every subset is both open and closed.
# The subsets of X are: {}, {1}, {2}, {1, 2}.
# In the discrete topology, all of these are closed sets.
all_subsets = [set(), {1}, {2}, {1, 2}]
print(f"\nWe equip X with the discrete topology. The closed sets are:")
for s in all_subsets:
    print(f"- {s or '{}'}")


# A "proper" closed subset is a closed subset Z such that Z is not equal to X.
proper_closed_sets = [s for s in all_subsets if s != X]
print(f"\nThe proper closed subsets of X are:")
for s in proper_closed_sets:
    print(f"- {s or '{}'}")


# Now, let's see if we can find two proper closed subsets, Z1 and Z2,
# whose union is X.
Z1 = {1}
Z2 = {2}
print(f"\nWe can choose two proper closed subsets: Z1 = {Z1} and Z2 = {Z2}.")

# Calculate the union
union_of_Z = Z1.union(Z2)

# We check if X = Z1 U Z2
print(f"The union of these two sets is {Z1} U {Z2} = {union_of_Z}.")

if union_of_Z == X:
    print(f"This union is equal to the original space X.")
    print("Since X is a union of two of its proper closed subsets, this 2-point space is NOT irreducible.")
else:
    print("The union is not equal to X.")

print("\nSummary of reasoning:")
print(" - For n=0 and n=1, all possible spaces are irreducible.")
print(" - For n=2, we have constructed a space that is not irreducible.")
print("Therefore, the smallest non-negative integer n for which a non-irreducible space exists is 2.")

# Final equation with each number printed.
print("\nThe final equation demonstrating this is:")
# The format "{1} U {2} = {1, 2}" shows each number.
z1_str = f"{{{list(Z1)[0]}}}"
z2_str = f"{{{list(Z2)[0]}}}"
x_list_sorted = sorted(list(X))
x_str = f"{{{x_list_sorted[0]}, {x_list_sorted[1]}}}"
print(f"{z1_str} U {z2_str} = {x_str}")
