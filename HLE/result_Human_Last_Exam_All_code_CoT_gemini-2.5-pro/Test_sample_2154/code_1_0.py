import math

# This script calculates the minimal order u_r(n) of the Picard-Fuchs
# differential equation for the Hamiltonian V(q) = 1/2 * (q^2 - q^n).
# The problem asks for the sequence of these orders for n from 3 to 12.

# The minimal order u_r(n) is given by 2*g, where g is the genus of the
# associated hyperelliptic curve y^2 = q^n - q^2 + 2*alpha.
# The genus is g = floor((n-1)/2).
# Thus, the formula for the order is u_r(n) = 2 * floor((n-1)/2).

# We calculate this for each n from 3 to 12.
results = []
n_range = range(3, 13)

for n in n_range:
    # The formula is 2 * floor((n-1)/2).
    # In Python, integer division `//` on positive numbers is equivalent
    # to taking the floor of the division.
    order = 2 * ((n - 1) // 2)
    results.append(order)

# The problem asks for the sequence of values {u_r(3), u_r(4), ..., u_r(12)}.
# To show the details of the final result, we first print each individual value
# that makes up the sequence.
print("The sequence of minimal orders {u_r(3), u_r(4), ..., u_r(12)} is calculated as follows:")
for n, order in zip(n_range, results):
    print(f"u_r({n}) = {order}")

# Finally, we print the complete sequence as a list, which is the answer.
print("\nThe complete sequence of values is:")
print(results)