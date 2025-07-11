import math

def S(m):
    """
    Calculates the sum of coefficients d_j for j from 0 to m.
    d_j = k+1 if j = 4k, and 0 otherwise.
    The sum is S(m) = sum_{k=0 to floor(m/4)} (k+1).
    """
    if m < 0:
        return 0
    # N is the upper limit for k, which is floor(m/4)
    N = m // 4
    # The sum of integers from 1 to N+1 is (N+1)(N+2)/2
    return (N + 1) * (N + 2) // 2

# Calculate the sum of coefficients up to the required degrees
s100 = S(100)
s98 = S(98)
s97 = S(97)
s95 = S(95)
s94 = S(94)

# The total rank is the sum of these values
total_rank = s100 + s98 + s97 + s95 + s94

# Print the equation with the calculated numbers, as requested
print("The total rank is the sum of coefficients of the Poincaré series up to degree 100.")
print("This can be calculated as: Total Rank = S(100) + S(98) + S(97) + S(95) + S(94)")
print(f"S(100) = {s100}")
print(f"S(98) = {s98}")
print(f"S(97) = {s97}")
print(f"S(95) = {s95}")
print(f"S(94) = {s94}")
print("\nThe final equation is:")
print(f"{s100} + {s98} + {s97} + {s95} + {s94} = {total_rank}")
print("\nThe total rank of A in degree <= 100 is:")
print(total_rank)