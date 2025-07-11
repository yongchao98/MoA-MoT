import sympy

# Step 1: Define the symbolic variable for the dimensionless conductance 'g'.
g = sympy.Symbol('g')

# Step 2: Define the universal probability distribution P(g) for a disordered
# Majorana wire at the critical point (Class D).
# P(g) = 1 / (2 * sqrt(g)) for g in the interval [0, 1].
P_g = 1 / (2 * sympy.sqrt(g))

# Step 3: Calculate the required statistical moments using symbolic integration.
# The n-th moment is given by the integral of g^n * P(g) from 0 to 1.

# Calculate the fourth moment, <g^4>, by setting n=4.
fourth_moment = sympy.integrate(g**4 * P_g, (g, 0, 1))

# Calculate the average value, <g>, which is the first moment (n=1).
average_value = sympy.integrate(g**1 * P_g, (g, 0, 1))

# Step 4: Compute the ratio of the fourth moment to the average value.
ratio = fourth_moment / average_value

# Step 5: Print the final equation, showing each calculated number as requested.
print("The fourth statistical moment is <g^4> = {}".format(fourth_moment))
print("The average value (first moment) is <g> = {}".format(average_value))
print("The ratio between the fourth moment and the average value is:")
print("<g^4> / <g> = ({}) / ({}) = {}".format(fourth_moment, average_value, ratio))
