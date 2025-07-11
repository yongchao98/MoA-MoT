import math

# This script calculates the critical time 'c' for the emergence of the
# giant component in the described variant of the Erdos-Renyi random graph model.

# Step 1: Determine the expected number of vertices N(t) at time t.
# Vertices appear according to a Poisson process with rate n.
# For a large n limit, the number of vertices N(t) is approximately n*t.
# N(t) ≈ n * t

# Step 2: Determine the expected number of edges M(t) at time t.
# The rate of edge formation at any time s is (N(s) choose 2) / n.
# For large n, N(s) ≈ n*s, so the rate is approximately ((n*s)^2 / 2) / n = (n * s^2) / 2.
# To find the total number of edges M(t), we integrate this rate from 0 to t:
# M(t) = Integral from 0 to t of (n * s^2 / 2) ds = n * t^3 / 6.

# Step 3: Calculate the average degree k(t) of the graph.
# The average degree is given by k(t) = 2 * M(t) / N(t).
# Substituting our expressions for N(t) and M(t):
# k(t) ≈ 2 * (n * t^3 / 6) / (n * t)
# k(t) ≈ t^2 / 3

# Step 4: Find the critical time 'c' by setting the average degree to 1.
# The giant component emerges when the average degree k(c) equals the threshold of 1.
# This gives us the equation to find c.

# The final equation is: c^2 / 3 = 1
# We can represent the components of this equation with variables:
exponent = 2
denominator = 3
threshold = 1

# Now, we solve this equation for 'c'.
# c^2 = 3 * 1
# c = sqrt(3)
c_squared = denominator * threshold
c_value = math.sqrt(c_squared)

print("The average degree k(t) as a function of time t is: k(t) = t^2 / 3")
print("The giant component emerges when the average degree k(c) reaches the critical threshold of 1.")
print("This leads to the final equation for the critical time c:")
print(f"c^{exponent} / {denominator} = {threshold}")

print("\nSolving the equation for c:")
print(f"c^{exponent} = {c_squared}")
print(f"c = sqrt({c_squared})")

print(f"\nThe exact value of c is the square root of 3.")
print(f"The numerical value of c is approximately: {c_value}")
