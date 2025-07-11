import math

# Time at which we want to find the position
t = 2 * math.sqrt(3)

# The trajectory x(t) is given by the sum of the trajectory in the accelerating frame, y(t),
# and the position of the frame itself, t^2/4.
# Based on the analysis, a key feature of the specified trajectory is that it reaches the origin
# of the accelerating frame at the specific time t = 2*sqrt(3).
# So, y(t) = 0 at t = 2*sqrt(3).

y_t = 0
x_t = y_t + (t**2) / 4

# The final equation is x(2*sqrt(3)) = 0 + (2*sqrt(3))^2 / 4
print(f"The position of the particle x(t) at t = 2*sqrt(3) is calculated by the equation:")
print(f"x({t:.4f}) = y({t:.4f}) + ({t:.4f}^2) / 4")
print(f"x({t:.4f}) = {y_t} + {t**2:.4f} / 4")
print(f"x({t:.4f}) = {x_t}")

# Final Answer
# print(f"<<<{x_t}>>>")