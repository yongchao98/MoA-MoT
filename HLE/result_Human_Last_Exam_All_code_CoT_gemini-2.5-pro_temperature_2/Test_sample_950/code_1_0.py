import numpy as np

# Define parameters for the Grassmannian Gr(k, n)
k = 3
n = 5

# --- Step 1: State the direct answer based on a mathematical theorem ---
print(f"The integral cohomology ring of the Grassmannian Gr({k}, {n}), denoted H*(Gr({k}, {n}); Z), is known to be torsion-free.")
print("This means its torsion subgroup contains only the zero element, i.e., it is the trivial group {0}.")
print("The rank of the trivial group {0} is 0.")
print("---")

# --- Step 2: Perform a related calculation for demonstration ---
# We calculate the Betti numbers of Gr(k, n), which are the ranks of the
# individual cohomology groups H^i(Gr(k,n); Z).

# By the duality property of Grassmannians, Gr(k, n) is diffeomorphic to Gr(n-k, n).
# The formula for the Poincare polynomial is simpler with the smaller index, k' = min(k, n-k).
k_calc = min(k, n - k)

# The Poincare polynomial P(t) for Gr(k, n) is given by the formula:
# P(t) = Product_{i=1 to k'} [ (1 - t^(n-i+1)) / (1 - t^i) ]
# For Gr(2, 5), this is [(1-t^5)*(1-t^4)] / [(1-t^1)*(1-t^2)].

# Calculate the numerator polynomial: Product_{i=1 to k'} (1 - t^(n-i+1))
p_num = np.poly1d([1])
for i in range(1, k_calc + 1):
    power = n - i + 1
    # Creates a polynomial like 1 - t^power
    poly_coeffs = [1] + [0] * (power - 1) + [-1]
    p_num = np.polymul(p_num, np.poly1d(poly_coeffs))

# Calculate the denominator polynomial: Product_{i=1 to k'} (1 - t^i)
p_den = np.poly1d([1])
for i in range(1, k_calc + 1):
    power = i
    poly_coeffs = [1] + [0] * (power - 1) + [-1]
    p_den = np.polymul(p_den, np.poly1d(poly_coeffs))

# The Poincare polynomial is the quotient of the two.
poincare_poly, remainder = np.polydiv(p_num, p_den)
# The coefficients are the Betti numbers. We round them to the nearest integer
# and reverse the list to have them in increasing order of index (b_0, b_1, ...).
betti_numbers = np.round(poincare_poly.coeffs).astype(int)[::-1]

print(f"For context, we can compute the Betti numbers (b_i) of Gr({k}, {n}), which represent the ranks of the (torsion-free) cohomology groups H^i.")
for i, b in enumerate(betti_numbers):
    print(f"b_{i} = {b}")

# The rank of the full cohomology ring is the sum of the Betti numbers.
# The "final equation" is presented here as the sum of Betti numbers.
total_rank = np.sum(betti_numbers)
equation_str = " + ".join(map(str, betti_numbers))
print(f"\nThe rank of the entire cohomology ring H*(Gr({k}, {n}); Z) is the sum:")
print(f"{equation_str} = {total_rank}")

# --- Step 3: Conclude with the final answer to the user's question ---
print("\nTherefore, the rank of the torsion subgroup is definitively:")
print("0")