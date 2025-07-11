import math

def calculate_u_r(n):
    """
    Calculates the minimal order of the Picard-Fuchs differential equation u_r(n)
    using the formula floor((n-1)/2).
    """
    return math.floor((n - 1) / 2)

# Calculate the sequence for n from 3 to 12
results = []
print("Calculating the values of u_r(n) for n = 3 to 12:")
for n in range(3, 13):
    u_r_n = calculate_u_r(n)
    results.append(u_r_n)
    print(f"u_r({n}) = {u_r_n}")

# Print the final sequence as a list
print("\nThe required set of values is:")
print(results)