# Let X be a smooth quintic hypersurface in CP^3.
# We want to find the rank of the third homotopy group pi_3(X).

# Step 1: Define the parameters of the hypersurface.
# d is the degree of the hypersurface.
d = 5
# n is the dimension of the ambient complex projective space CP^n.
n = 3

# Step 2: Use topological properties to find the Betti numbers of X.
# b_k = rank(H_k(X)).
# By the Lefschetz Hyperplane Theorem, pi_1(X) = 0, so H_1(X) = 0 and b_1 = 0.
b1 = 0
# X is a compact orientable 4-manifold (complex dimension 2).
# By Poincare Duality, b_3 = b_1 = 0.
b3 = 0
# X is connected, so b_0 = 1.
b0 = 1
# By Poincare Duality, b_4 = b_0 = 1.
b4 = 1

# Step 3: Calculate the second Betti number, b_2, using the Euler characteristic.
# The Euler characteristic chi(X) is given by the sum of Betti numbers:
# chi(X) = b_0 - b_1 + b_2 - b_3 + b_4 = 1 - 0 + b_2 - 0 + 1 = 2 + b_2.
# The Euler characteristic can also be computed from the degree d using the formula:
# chi(X) = d * (d^2 - (n+1)*d + (n+1)*n/2 + 1)
# For n=3, this simplifies to chi(X) = d * (d^2 - 4*d + 6).
chi = d * (d**2 - 4*d + 6)
# Now we can solve for b_2: b_2 = chi - 2.
b2 = chi - 2

print(f"The Euler characteristic of X is chi(X) = {d} * ({d**2} - 4*{d} + 6) = {chi}.")
print(f"The second Betti number is b_2 = chi(X) - 2 = {chi} - 2 = {b2}.")

# Step 4: Relate the Betti numbers to homotopy groups.
# Since pi_1(X) = 0, the Hurewicz theorem states pi_2(X) is isomorphic to H_2(X).
# The rank of pi_2(X) is therefore b_2.
rank_pi2 = b2

# Step 5: Calculate the rank of pi_3(X).
# For a simply connected space X, Whitehead's theorem states that if H_3(X) = 0 (which is true as b_3=0),
# then pi_3(X) is isomorphic to Gamma(pi_2(X)), where Gamma is Whitehead's universal quadratic functor.
# The rank of Gamma(Z^k) is given by the formula k * (k - 1) / 2.
k = rank_pi2
rank_pi3 = k * (k - 1) // 2

print(f"\nThe rank of pi_2(X) is {k}.")
print(f"The rank of pi_3(X) is given by the formula rank(Gamma(Z^k)) = k * (k - 1) / 2.")
print(f"Rank(pi_3(X)) = {k} * ({k} - 1) / 2 = {k} * {k-1} / 2 = {rank_pi3}")
