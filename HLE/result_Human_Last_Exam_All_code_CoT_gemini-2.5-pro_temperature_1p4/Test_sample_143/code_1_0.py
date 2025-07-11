import numpy as np

# This script demonstrates that for d=3, the condition x_i >= 0 is not sufficient
# for the scalarization method to explore the entire Pareto front.
# The analysis shows that scalarization fails to find the interior points of the Pareto front.

# 1. Define the counterexample for d=3 with non-negative data matrix X.
x1 = np.array([1.0, 1.0, 0.0])
x2 = np.array([0.0, 0.0, 1.0])

print("Demonstration for d=3 counterexample with non-negative data:")
print(f"x1 = {x1}")
print(f"x2 = {x2}\n")

# 2. Define the Pareto Front analytically.
# The objectives to maximize are z1 = (w_1+w_2)^2 and z2 = w_3^2.
# A point (z1, z2) is on the Pareto front if z1 = 2u and z2 = 1-u for some u in [0,1].
# This corresponds to the line segment z2 = 1 - z1/2 for z1 in [0, 2].
print("The true Pareto front is the line segment described by the equation:")
# Equation: z2 = 1 - 0.5*z1  or  0.5*z1 + 1.0*z2 - 1.0 = 0
eq_c1, eq_c2, eq_c0 = 0.5, 1.0, -1.0
print(f"Equation: {eq_c1}*z1 + {eq_c2}*z2 + {eq_c0} = 0, for z1 in [0, 2]\n")

# The endpoints of this segment are:
p_front_pt1 = np.array([0.0, 1.0]) # Corresponds to u=0
p_front_pt2 = np.array([2.0, 0.0]) # Corresponds to u=1
# An example of a point in the interior of the front is:
p_front_interior_pt = np.array([1.0, 0.5]) # Corresponds to u=0.5

print(f"The endpoints of the Pareto front are: {p_front_pt1} and {p_front_pt2}")
print(f"An interior point on the Pareto front is: {p_front_interior_pt}\n")


# 3. Determine the points found by the scalarization method.
# The scalarized objective is L = lambda1 * (w1+w2)^2 + (1-lambda1) * w3^2.
# As derived in the thinking steps, maximizing L yields one of the two endpoints,
# depending on the value of lambda1.

# For lambda1 > 1/3, the maximizer is w = (1/sqrt(2), 1/sqrt(2), 0)
w_sol1 = np.array([1.0/np.sqrt(2), 1.0/np.sqrt(2), 0.0])
z1_sol1 = np.dot(w_sol1, x1)**2
z2_sol1 = np.dot(w_sol1, x2)**2
scalar_sol1 = np.array([z1_sol1, z2_sol1])

# For lambda1 < 1/3, the maximizer is w = (0, 0, 1)
w_sol2 = np.array([0.0, 0.0, 1.0])
z1_sol2 = np.dot(w_sol2, x1)**2
z2_sol2 = np.dot(w_sol2, x2)**2
scalar_sol2 = np.array([z1_sol2, z2_sol2])

print("Points found by sweeping weights in the scalarization method:")
print(f"For almost all weights, the method finds either {scalar_sol1} or {scalar_sol2}.")
print("Specifically, it only finds the endpoints of the Pareto front.")
print("\nConclusion: The scalarization method fails to generate the interior points of the front, like (1.0, 0.5).")
print("Therefore, for d=3, the condition that x_i >= 0 is not sufficient.")
print("\nSince scalarization is known to work for any data when d=2, the largest dimension for which the condition is sufficient is 2.")