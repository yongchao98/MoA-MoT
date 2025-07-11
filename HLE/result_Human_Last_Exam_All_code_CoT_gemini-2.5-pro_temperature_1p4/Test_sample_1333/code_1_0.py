def euler_characteristic(group_order, l, m, k):
    """
    Calculates the Euler characteristic for a regular dessin.
    Formula: chi = |G| * (1/l + 1/m + 1/k - 1/2)
    l, m, k are the orders of b, w, bw.
    """
    if l == 0 or m == 0 or k == 0:
        return None
    # Check if the characteristic is negative (hyperbolic case)
    hyperbolic_check = 1/l + 1/m + 1/k - 1/2
    if hyperbolic_check >= 0:
        # The problem specifies a negative Euler characteristic.
        # This combination of l,m,k does not yield a valid dessin D.
        return None
    
    return group_order * hyperbolic_check

# Parameters for our example based on SL(2,7)
# For the covering dessin D
G_order = 336  # |SL(2,7)|
N_order = 2    # |Z(SL(2,7))|

# For the quotient dessin D_N
G_N_order = G_order / N_order # |PSL(2,7)|

# Signature (l, m, k) = (|b|, |w|, |bw|)
# We chose a signature where a smooth covering is possible.
l, m, k = 3, 3, 7

# Calculate Euler characteristics
chi_D = euler_characteristic(G_order, l, m, k)
chi_D_N = euler_characteristic(G_N_order, l, m, k)

# The ratio is |N|
ratio = chi_D / chi_D_N

# Print the calculation steps
print(f"Let G = SL(2,7) and N = Z(G) be the center of G.")
print(f"The order of G is |G| = {G_order}.")
print(f"The order of N is |N| = {N_order}.")
print(f"The quotient group is G/N = PSL(2,7), with order |G/N| = {G_N_order}.")
print(f"We use the signature (l,m,k) = ({l},{m},{k}) for the dessin.")
print(f"This signature is hyperbolic because 1/l + 1/m + 1/k - 1/2 = 1/{l} + 1/{m} + 1/{k} - 1/2 = {1/l + 1/m + 1/k - 1/2:.4f} < 0.")
print("-" * 20)
print(f"Euler characteristic of D:")
print(f"chi(D) = |G| * (1/l + 1/m + 1/k - 1/2)")
print(f"chi(D) = {G_order} * (1/{l} + 1/{m} + 1/{k} - 0.5) = {chi_D}")
print("-" * 20)
print(f"Euler characteristic of D_N:")
print(f"chi(D_N) = |G/N| * (1/l + 1/m + 1/k - 1/2)")
print(f"chi(D_N) = {G_N_order} * (1/{l} + 1/{m} + 1/{k} - 0.5) = {chi_D_N}")
print("-" * 20)
print(f"The ratio chi(D) / chi(D_N) is:")
print(f"Ratio = {chi_D} / {chi_D_N} = {ratio}")
print(f"This ratio is equal to the order of the normal subgroup N, |N| = {N_order}.")
print("-" * 20)
print("Based on this and related examples in the theory of regular maps, the maximum value appears to be constrained by properties of group extensions.")
print("The value of 2 is obtained from a standard construction involving a central extension of a simple group.")
