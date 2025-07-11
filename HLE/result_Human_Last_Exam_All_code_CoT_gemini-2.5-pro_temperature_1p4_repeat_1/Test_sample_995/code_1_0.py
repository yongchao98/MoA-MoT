import math

# The problem asks us to find the minimized probability P(A wins)
# and then calculate floor(1 / P(A wins)).

# From the derivation, the minimum probability is given by the formula:
# P_min = (2/pi) * arctan(K) - 1/2
# where K = 2 + pi/4

pi_val = math.pi
K = 2 + pi_val / 4
p_min = (2 / pi_val) * math.atan(K) - 0.5
result = math.floor(1 / p_min)

# Output the equation and the values of its components as requested.
print(f"The equation for the minimized probability is P(A wins) = (2/pi) * arctan(2 + pi/4) - 1/2")
print(f"The value of pi is: {pi_val}")
print(f"The value of K (2 + pi/4) is: {K}")
print(f"The minimized probability P(A wins) is: {p_min}")
print(f"The reciprocal 1/P(A wins) is: {1/p_min}")
print(f"The floor of the reciprocal is:")
print(result)
