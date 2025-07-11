import math

# This script calculates the constant 'c', the time of emergence of the giant
# connected component in the specified random graph model.

# The solution follows these steps:
# 1. Model the number of vertices V(t) and edges E(t) over time.
#    - Vertex arrival rate is n, so in the n->inf limit, V(t) = n*t.
#    - Edge arrival rate is 1/n for each possible pair. The number of pairs
#      at time t is V(t)-choose-2, which is approx. (n*t)^2 / 2.
#    - The rate of edge growth is dE/dt = (1/n) * (V(t) choose 2)
#      dE/dt = (1/n) * (n*t)^2 / 2 = (n * t^2) / 2.

# 2. Integrate dE/dt to find E(t).
#    E(t) = integral from 0 to t of (n * s^2 / 2) ds
#    E(t) = (n / 2) * [s^3 / 3] from 0 to t = n * t^3 / 6.

# 3. Calculate the average degree k(t) = 2 * E(t) / V(t).
#    k(t) = 2 * (n * t^3 / 6) / (n * t)
#    k(t) = t^2 / 3.

# 4. The giant component emerges when the average degree k(c) = 1.
#    This gives the equation we need to solve: c^2 / 3 = 1.

# --- Calculation ---
# Here we solve the equation for c, printing each step.

# The equation for the critical time c is:
eq_lhs_numerator_var = "c^2"
eq_denominator = 3
eq_rhs = 1
print("The average degree k(t) of the graph at time t is t^2 / 3.")
print(f"The giant component emerges when the average degree is {eq_rhs}.")
print("This gives us the equation:")
print(f"{eq_lhs_numerator_var} / {eq_denominator} = {eq_rhs}")

# Solve for c^2
c_squared = eq_rhs * eq_denominator
print("\nSolving for c^2:")
print(f"c^2 = {eq_rhs} * {eq_denominator}")
print(f"c^2 = {c_squared}")

# Solve for c
c = math.sqrt(c_squared)
print("\nSolving for c:")
print(f"c = sqrt({c_squared})")
print(f"\nThe exact value of c is the square root of {c_squared}.")
print(f"The numerical value of c is approximately: {c}")

<<<sqrt(3)>>>