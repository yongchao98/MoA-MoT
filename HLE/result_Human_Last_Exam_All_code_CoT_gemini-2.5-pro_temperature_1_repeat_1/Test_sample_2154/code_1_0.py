# This script calculates the minimal order of the Picard-Fuchs differential equation u_r(n)
# for the Hamiltonian V(q) = 1/2 * (q^2 - q^n) for n from 3 to 12.

# The formula for u_r(n) depends on the parity of n:
# - If n is odd, u_r(n) = n - 1.
# - If n is even, u_r(n) = n / 2, due to the symmetry of the potential V(q).

# We will generate the list of values {u_r(3), u_r(4), ..., u_r(12)}.
u_r_values = []
for n in range(3, 13):
    if n % 2 == 0:  # n is even
        # For even n, the order is n/2.
        order = n // 2
    else:  # n is odd
        # For odd n, the order is n-1.
        order = n - 1
    u_r_values.append(order)

# The problem asks to output each number. The clearest way is to print the list.
print(u_r_values)