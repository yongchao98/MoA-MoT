import math

def u_r(n):
    """
    Calculates the minimal order of the Picard-Fuchs differential equation
    for the Hamiltonian V(q) = 1/2 * (q^2 - q^n).
    
    The formula is derived from a known result in mathematics, which simplifies to n-2.
    u_r(n) = floor((n-1)/2) + floor((n-2)/2) = n-2
    """
    if n < 3:
        raise ValueError("The formula is defined for n >= 3.")
    
    # The simplified formula is u_r(n) = n - 2
    return n - 2

# The range of n is from 3 to 12, inclusive.
start_n = 3
end_n = 12

# We will store the results in a list to display at the end.
results = []

# Calculate and print u_r(n) for each n in the specified range.
print(f"Calculating u_r(n) for n from {start_n} to {end_n}:")
for n in range(start_n, end_n + 1):
    order = u_r(n)
    results.append(order)
    print(f"u_r({n}) = {n} - 2 = {order}")

# Final answer as a list of values for {u_r(3), u_r(4), ..., u_r(12)}
print("\nThe set of values is:")
print(results)