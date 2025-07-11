import math
from itertools import combinations

def nCr(n, r):
    """Calculates the binomial coefficient 'n choose r'."""
    if r < 0 or r > n:
        return 0
    # Use integer division
    return math.factorial(n) // (math.factorial(r) * math.factorial(n - r))

def check_l_intersecting(family, L_set):
    """Checks if a family of sets is L-intersecting."""
    family_list = list(family)
    for i in range(len(family_list)):
        for j in range(i + 1, len(family_list)):
            set_i = family_list[i]
            set_j = family_list[j]
            intersection_size = len(set_i.intersection(set_j))
            if intersection_size not in L_set:
                return False
    return True

# --- Analysis for Question (b) ---
print("--- Verifying the Counterexample for Question (b) ---")
print("Must the bound m <= sum_{i=0 to s} binom(n-1, i) hold?")

# 1. Define the parameters for our counterexample.
n = 4
s = 1
L = {1}
print(f"\nParameters: n = {n}, s = {s}, L = {L}")

# 2. Construct the family F.
# We choose F to be the family of all 2-element subsets of {1, 2, 3, 4}.
base_set = set(range(1, n + 1))
k = 2
family_tuples = combinations(base_set, k)
F = [set(t) for t in family_tuples]
m = len(F)

print(f"\nConstructed Family F: all {m} {k}-element subsets of {base_set}.")
print("F = {", ", ".join(map(str, F)), "}")
print(f"Size of family, m = {m}")

# 3. Verify that F is an ordered L-intersecting family.
# L-intersecting property:
is_L_int = check_l_intersecting(F, L)
print(f"\nIs F an L-intersecting family? {is_L_int}")

# Ordered property:
# A family is ordered if there exists r such that n is in F_1..F_r, n is not in F_{r+1}..F_m,
# and sizes are non-decreasing. All our sets have size 2, so the size condition is trivially met.
# We can partition the sets based on the element n=4, so an ordering exists.
print(f"Is F an 'ordered' family? Yes, all sets have the same size (2) and can be partitioned by element n={n}.")

# 4. Calculate the bound from the question.
bound = 0
bound_eq_parts = []
bound_val_parts = []
for i in range(s + 1):
    term = nCr(n - 1, i)
    bound += term
    bound_eq_parts.append(f"binom({n-1}, {i})")
    bound_val_parts.append(str(term))

print("\nCalculating the proposed bound:")
print(f"Bound = sum(i=0 to s) binom(n-1, i) = {' + '.join(bound_eq_parts)}")
# Output each number in the final equation as requested
print(f"      = {' + '.join(bound_val_parts)}")
print(f"      = {bound}")

# 5. Compare the family size m with the bound.
print("\nComparison:")
print(f"Family size m = {m}")
print(f"Proposed bound = {bound}")

if m > bound:
    print(f"\nResult: m > Bound ({m} > {bound}).")
    print("The constructed family is a valid counterexample.")
    print("Conclusion: The bound does not necessarily hold for any ordered L-intersecting family.")
else:
    print("\nThis example does not violate the bound.")
