import math

# The emergence of the giant component in a random graph occurs when the
# average degree of the graph reaches 1. Let 'c' be the time this happens.
# The average degree is defined as 2 * E[Edges] / E[Vertices], where E denotes
# the expected value. Our goal is to find expressions for the expected number
# of vertices and edges at time 'c'.

# --- Step 1: Expected number of vertices at time c ---
# Vertices arrive via a Poisson process with rate 'n'.
# The expected number of vertices at time t is E[V(t)] = n*t.
# So, at time c, E[V(c)] = n*c.
n_vertices_str = "n * c"
print(f"Let c be the critical time for the emergence of the giant component.")
print(f"Let n be the rate of vertex arrival.")
print(f"The expected number of vertices at time c is E[V(c)] = {n_vertices_str}")
print("-" * 60)

# --- Step 2: Expected number of edges at time c ---
# At any given time t, the number of vertices, k, follows a Poisson distribution
# with mean (lambda) = n*t.
# The number of possible edges is k_choose_2 = k*(k-1)/2.
# For a Poisson distribution, the expected value of k*(k-1) is lambda^2.
# So, the expected number of vertex pairs at time t is E[k*(k-1)/2] = (n*t)^2 / 2.
# Edges appear between any existing pair at a rate of 1/n.
# Therefore, the rate of edge formation at time t is:
# Rate(t) = (1/n) * E[vertex pairs] = (1/n) * (n*t)^2 / 2 = n * t^2 / 2.
# To find the expected number of edges at time c, we integrate this rate from 0 to c.
# E[E(c)] = integral from 0 to c of (n * t^2 / 2) dt = n * c^3 / 6.
n_edges_str = "n * c**3 / 6"
print("The expected rate of edge formation at time t is n * t**2 / 2.")
print("Integrating this rate from 0 to c gives the expected number of edges at time c:")
print(f"E[E(c)] = {n_edges_str}")
print("-" * 60)

# --- Step 3: Set average degree to 1 and solve for c ---
# Average Degree = 2 * E[E(c)] / E[V(c)]
# The giant component emerges when the Average Degree = 1.
# 1 = 2 * (n * c^3 / 6) / (n * c)
# 1 = (n * c^3 / 3) / (n * c)
# 1 = c^2 / 3
avg_degree_str = "c**2 / 3"
critical_value = 1
print("The giant component emerges when the average degree is 1.")
print("Average Degree = 2 * E[E(c)] / E[V(c)] = c**2 / 3")
print(f"Setting the average degree to 1, we get the equation:")
print(f"{avg_degree_str} = {critical_value}")
print("")

# Solving the equation c^2 / 3 = 1
# c^2 = 3
c_squared_val = 3
print("Solving for c:")
print(f"c**2 = {c_squared_val}")

# c = sqrt(3)
c = math.sqrt(c_squared_val)
print(f"c = sqrt({c_squared_val})")
print("")
print(f"The exact value of c is sqrt(3), which is approximately {c}")
