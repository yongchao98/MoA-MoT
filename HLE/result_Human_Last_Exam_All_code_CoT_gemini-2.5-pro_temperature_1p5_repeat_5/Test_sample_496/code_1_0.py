import math

def S(N):
    """
    Calculates the sum of ranks of H^*(BSO(4)) up to degree N.
    The Poincare series of H^*(BSO(4)) is 1/((1-t^4)^2), which expands to
    sum_{j=0 to infinity} (j+1)*t^(4j).
    S(N) computes sum_{k=0 to N} [t^k]P(t), which is sum_{j=0 to floor(N/4)} (j+1).
    """
    M = N // 4
    # The sum is 1 + 2 + ... + (M+1)
    return (M + 1) * (M + 2) // 2

# Calculate the components of the total sum
s100 = S(100)
s97 = S(97)
s94 = S(94)

# The total rank is given by the sum of coefficients of the Poincare series
# (1 + 2*t^3 + t^6) * (sum_{j=0 to infinity} (j+1)*t^(4j))
# up to degree 100.
total_rank = s100 + 2 * s97 + s94

print(f"The total rank is a sum of three components derived from the Poincare series.")
print(f"Component 1 (from the '1' term): Sum of ranks of H^*(BSO(4)) up to deg 100 is S(100) = {s100}")
print(f"Component 2 (from the '2t^3' term): 2 * Sum of ranks of H^*(BSO(4)) up to deg 97 is 2 * S(97) = 2 * {s97} = {2*s97}")
print(f"Component 3 (from the 't^6' term): Sum of ranks of H^*(BSO(4)) up to deg 94 is S(94) = {s94}")
print(f"\nThe equation for the total rank is: {s100} + 2 * {s97} + {s94}")
print(f"Total rank = {total_rank}")
