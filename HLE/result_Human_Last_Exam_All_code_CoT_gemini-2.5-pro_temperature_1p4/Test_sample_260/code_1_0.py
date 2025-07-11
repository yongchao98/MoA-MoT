from fractions import Fraction

# Step 1 & 2: Define the orders of the fundamental groups.
orders = [5, 8, 2]
num_groups = len(orders)

# Step 4: Calculate the Euler characteristics.

# Calculate chi(G) for G = Z_5 * Z_8 * Z_2.
# The formula is chi(G) = sum(chi(Gi)) - (n-1).
# For a finite cyclic group Z_m, chi(Z_m) = 1/m.
chi_Gi_list = [Fraction(1, n) for n in orders]
chi_G = sum(chi_Gi_list) - (num_groups - 1)

# Calculate chi(Q) for Q = Z_5 + Z_8 + Z_2.
# For a finite group Q, chi(Q) = 1/|Q|.
order_Q = 1
for n in orders:
    order_Q *= n
chi_Q = Fraction(1, order_Q)

# The Euler characteristic of the kernel K (a free group of rank r) is chi(K) = 1 - r.
# We use the formula chi(G) = chi(K) * chi(Q), which is chi(G) = (1-r) * chi(Q).

# Step 5: Solve for the rank r.
# r = 1 - chi(G) / chi(Q)
chi_K = chi_G / chi_Q
rank_r = 1 - chi_K

# Output the explanation and the final equation.
print("The problem reduces to solving for the rank 'r' in the group-theoretic equation involving Euler characteristics:")
print("chi(G) = chi(K) * chi(Q)")
print("where G = pi_1(Y), K is the kernel (a free group of rank r), and Q = H_1(Y).")
print("\nWe calculate the components:")
print(f"chi(G) = (1/5 + 1/8 + 1/2) - (3 - 1) = {chi_G}")
print(f"chi(Q) = 1 / (5 * 8 * 2) = {chi_Q}")
print(f"chi(K) = 1 - r")
print("\nSubstituting these values into the formula gives the equation:")
print(f"{chi_G} = (1 - r) * {chi_Q}")
print("\nSolving for the rank r:")
print(f"1 - r = {chi_G} / {chi_Q} = {chi_K}")
print(f"r = 1 - ({chi_K}) = {rank_r}")

print(f"\nThe rank of the kernel as a free group is {int(rank_r)}.")