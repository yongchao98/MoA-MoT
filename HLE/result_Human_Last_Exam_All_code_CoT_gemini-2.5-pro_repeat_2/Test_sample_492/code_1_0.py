import math

# We want to find the critical time 'c' for the emergence of the giant component.

# Step 1: Expected number of vertices at time t
# Vertices appear at a rate of n. This is a Poisson process.
# At a given time t, the expected number of vertices V(t) is:
# E[V(t)] = n * t
# For large n, the actual number of vertices is sharply concentrated around its mean.
# So, we can approximate V(t) ≈ n * t.
print("Step 1: Determine the expected number of vertices.")
print("Let 'c' be the critical time we are looking for.")
print("At time c, the expected number of vertices is V(c) ≈ n * c.\n")

# Step 2: Average edge probability at time t
# Consider any two vertices, u and v, present at time c. Their arrival times, T_u and T_v,
# are independent and uniformly distributed in the interval [0, c].
# An edge can form between them at rate 1/n only after both have arrived, i.e., in the time
# window [max(T_u, T_v), c]. The duration of this window is c - max(T_u, T_v).
# The probability of an edge forming in this window is p_uv = 1 - exp(-(c - max(T_u, T_v))/n).
# For large n, the term in the exponent is small, so we can use the approximation 1 - e^(-x) ≈ x.
# p_uv ≈ (c - max(T_u, T_v)) / n
# To get the average probability 'p' for any pair, we take the expectation over T_u and T_v.
# p = E[p_uv] ≈ E[(c - max(T_u, T_v)) / n] = (c - E[max(T_u, T_v)]) / n
# The expectation E[max(T_u, T_v)] for two i.i.d. uniform variables on [0, c] is 2*c/3.
# Let's show this calculation:
# For Z = max(T_u, T_v), the CDF is F_Z(z) = P(Z <= z) = P(T_u <= z)P(T_v <= z) = (z/c)^2.
# The PDF is f_Z(z) = d/dz F_Z(z) = 2*z/c^2.
# E[Z] = ∫[0 to c] z * (2*z/c^2) dz = (2/c^2) * [z^3/3]_0^c = (2/c^2)*(c^3/3) = 2*c/3.
e_max_t = "2*c/3"
# So, p ≈ (c - 2*c/3) / n = (c/3) / n.
print("Step 2: Calculate the average edge probability 'p'.")
print(f"The average probability 'p' between any two vertices at time 'c' is approximately (c - E[max(T_u,T_v)])/n.")
print(f"The expectation of the maximum of two uniform arrival times in [0,c] is E[max(T_u,T_v)] = {e_max_t}.")
print("This gives p ≈ (c - 2*c/3) / n = (c/3) / n.\n")

# Step 3: Average degree at time t
# The graph at time c resembles a G(N, p) random graph with N=V(c) vertices.
# The average degree λ is approximately N * p.
# λ ≈ (n * c) * ((c/3) / n) = c^2 / 3
print("Step 3: Calculate the average degree λ.")
print("The average degree λ ≈ V(c) * p ≈ (n * c) * ((c/3) / n).")
print("λ ≈ c^2 / 3.\n")

# Step 4: Apply the critical condition
# The giant component emerges when the average degree λ equals 1.
# We set up the equation to solve for c.
numerator = 1
denominator = 3
rhs = 1
print("Step 4: Find 'c' by setting the average degree to 1.")
print("The critical condition for the emergence of a giant component is λ = 1.")
print("The final equation is:")
print(f"c^2 / {denominator} = {rhs}")

# Solving the equation
c_squared = denominator * rhs
c_numeric = math.sqrt(c_squared)
print("\nSolving for c:")
print(f"c^2 = {c_squared}")
print(f"c = sqrt({c_squared})")
print(f"The exact value of c is sqrt(3), which is approximately {c_numeric}.")
