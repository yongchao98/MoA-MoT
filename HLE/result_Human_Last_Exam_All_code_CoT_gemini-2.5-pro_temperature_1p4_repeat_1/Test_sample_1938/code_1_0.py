import math

# We want to find the minimum q (q_0) such that the probability of Alice winning (p_A) is at least 0.95.
# Based on the analysis of the game, p_A and q are related by the following equation:
# q^4 * p_A^3 = 1 - (1 - p_A)^(1/3)
# We set p_A = 0.95 to find the corresponding q_0.

p_A = 0.95

# The final equation to solve for q_0 is: q_0^4 * (0.95)^3 = 1 - (1-0.95)^(1/3)
# Let's compute each part of the equation.

# Right Hand Side (RHS) of the equation
rhs = 1 - (1 - p_A)**(1/3)

# The p_A term on the Left Hand Side (LHS)
p_A_cubed = p_A**3

# Now, solve for q_0^4
# q_0^4 = RHS / p_A_cubed
q_pow_4 = rhs / p_A_cubed

# Finally, q_0 is the 4th root
q0 = q_pow_4**(1/4)

# The problem asks for the value of floor(100 * q_0)
result = math.floor(100 * q0)

# Print out the calculation steps with the numbers filled in
print("We are solving the equation: q_0^4 * p_A^3 = 1 - (1 - p_A)^(1/3)")
print(f"Setting p_A = {p_A}:\n")
print(f"The right-hand side is 1 - (1 - {p_A})^(1/3):")
print(f"1 - ({1-p_A})^(1/3) = {rhs}\n")
print(f"The p_A term on the left-hand side is ({p_A})^3:")
print(f"({p_A})^3 = {p_A_cubed}\n")
print("Solving for q_0^4:")
print(f"q_0^4 = {rhs} / {p_A_cubed} = {q_pow_4}\n")
print("Solving for q_0:")
print(f"q_0 = ({q_pow_4})^(1/4) = {q0}\n")
print(f"The value of 100 * q_0 is: {100 * q0}\n")
print(f"The final answer, floor(100 * q_0), is:")
print(result)