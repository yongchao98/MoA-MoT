import sys

# Plan:
# 1. Use the isomorphism G(3,5) ~ G(2,5). Let k=2, n=5. Dimension d = k(n-k) = 6.
# 2. Find the mod 2 Betti numbers (beta_i) from Schubert cell counts.
# 3. Find the rational Betti numbers (b_i) from the rational cohomology ring structure.
# 4. Use the formula beta_i = b_i + t_i + t_{i-1} to compute the ranks of homology torsion (t_i).
# 5. The total rank of cohomology torsion is the sum of all t_i.

# For G(2,5), k=2, n=5, dimension d=6.

# Step 1: Mod 2 Betti numbers (beta_i)
# These are the number of partitions (lambda_1, lambda_2) with 3 >= lambda_1 >= lambda_2 >= 0
# of a given size i = lambda_1 + lambda_2.
# i=0: (0,0) -> 1
# i=1: (1,0) -> 1
# i=2: (2,0), (1,1) -> 2
# i=3: (3,0), (2,1) -> 2
# i=4: (3,1), (2,2) -> 2
# i=5: (3,2) -> 1
# i=6: (3,3) -> 1
beta = [1, 1, 2, 2, 2, 1, 1]
print("Step 1: The Betti numbers with Z/2Z coefficients (beta_i) are:")
print(f"beta = {beta}\n")

# Step 2: Rational Betti numbers (b_i)
# H*(G(2,5); Q) = Q[p_1]/(p_1^2) where deg(p_1) = 4.
# So, b_0=1, b_4=1, and all other b_i are 0.
# The manifold is non-orientable, so b_d = b_6 = 0.
b = [1, 0, 0, 0, 1, 0, 0]
print("Step 2: The rational Betti numbers (b_i) are:")
print(f"b = {b}\n")

# Step 3: Compute ranks of homology torsion subgroups (t_i)
dim = 6
t = [0] * (dim + 1)
print("Step 3: Calculating the ranks of the integral homology torsion subgroups (t_i) using beta_i = b_i + t_i + t_{i-1}:")
for i in range(dim + 1):
    t_prev = t[i-1] if i > 0 else 0
    t[i] = beta[i] - b[i] - t_prev
    if t[i] < 0:
        print(f"Error: Calculation resulted in a negative rank at i={i}. Aborting.")
        sys.exit(1)
    print(f"i={i}: t_{i} = beta_{i} - b_{i} - t_{{{i-1}}} = {beta[i]} - {b[i]} - {t_prev} = {t[i]}")
print(f"\nThe ranks of homology torsion are:\nt = {t}\n")

# Step 4: Compute the total rank of the torsion subgroup of the integral cohomology ring.
# The rank of the torsion part of H^i is t_{i-1}.
# The total rank is the sum of these ranks: sum(t_i for i=0 to d).
total_rank = sum(t)
print("Step 4: The rank of the torsion subgroup of the integral cohomology ring is the sum of the ranks of the torsion part of each cohomology group H^i.")
print("The rank of Torsion(H^i) is t_{i-1}. The total rank is the sum of all t_j for j=0 to 6.")
equation_parts = [str(val) for val in t]
equation_str = " + ".join(equation_parts)
print(f"The calculation is: {equation_str} = {total_rank}")

print(f"\nThe rank of the torsion subgroup of the integral cohomology ring of the space of 3-subspaces of R^5 is {total_rank}.")
