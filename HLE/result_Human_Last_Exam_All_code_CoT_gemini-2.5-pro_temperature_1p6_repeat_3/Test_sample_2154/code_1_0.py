import math

def calculate_ur(n):
    """
    Calculates the minimal order of the Picard-Fuchs differential equation
    u_r(n) using the formula 2 * floor((n-1)/2).
    """
    # Using integer division // is equivalent to math.floor() for positive results.
    return 2 * ((n - 1) // 2)

# We need to find the values for n from 3 to 12.
n_values = range(3, 13)

results = []
print("Calculating the sequence {u_r(3), u_r(4), ..., u_r(12)}:")
for n in n_values:
    ur_n = calculate_ur(n)
    results.append(ur_n)
    # The prompt asks to output each number in the final equation.
    # We interpret this as showing the result for each n.
    print(f"u_r({n}) = 2 * floor(({n}-1)/2) = {ur_n}")

print("\nThe complete sequence is:")
print(results)