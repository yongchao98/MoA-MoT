import math

# Problem parameters
n = 4  # Dimension of the vector space D over F_2
q = 2  # Size of the field k

# --- Step 1: Calculate the order of GL(n, q) ---

# Calculate each term in the product formula
terms = [q**n - q**i for i in range(n)]

# Calculate the full order of the group
order_gl = 1
for term in terms:
    order_gl *= term

# Display the calculation for the order of GL(n, q)
term_strings = [str(t) for t in terms]
print(f"The structure of the defect group D is a {n}-dimensional vector space over F_{q}.")
print(f"The automorphism group Aut(D) is isomorphic to GL({n}, {q}).")
print("The order of GL(n, q) is the product of (q^n - q^i) for i from 0 to n-1.")
print(f"|GL({n}, {q})| = {' * '.join(term_strings)}")
print(f"|GL({n}, {q})| = {order_gl}")
print("-" * 30)

# --- Step 2: Find the largest odd divisor ---

# The inertial quotient E must be a 2'-group (odd order), so its order must be
# the odd part of |GL(n, q)|. We find this by dividing out all factors of 2.
order_sylow_2 = 1
temp_order = order_gl
while temp_order % 2 == 0 and temp_order != 0:
    temp_order //= 2
    order_sylow_2 *= 2

highest_order_E = order_gl // order_sylow_2

# Display the final calculation
print("The order of the inertial quotient E must be odd.")
print("The highest possible order for E is the largest odd divisor of |GL(4, 2)|.")
print("This is calculated by dividing the total order by the order of its Sylow 2-subgroup.")
print(f"Final Calculation: {order_gl} / {order_sylow_2} = {highest_order_E}")

print("\nThus, the highest possible order for the inertial quotient E is:", highest_order_E)
