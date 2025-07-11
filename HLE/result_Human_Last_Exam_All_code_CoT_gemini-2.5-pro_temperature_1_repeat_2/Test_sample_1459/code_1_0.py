import math

# The Gromov-Hausdorff distance 'd' is the solution to cos(1/(2d)) + cos(pi/(2d)) = 0.
# The analytical solution is d = (pi + 1) / (2 * pi).

pi = math.pi
d = (pi + 1) / (2 * pi)

# The equation is derived from setting u = 1/(2d) in cos(u) + cos(pi*u) = 0
# which gives u(pi+1) = pi, so u = pi / (pi+1)
# 1/(2d) = pi / (pi+1)
# 2d = (pi+1)/pi
# d = (pi+1)/(2*pi)
# We can also write it as d = 1/2 + 1/(2*pi)

val_1_div_2 = 1/2
val_1_div_2pi = 1 / (2 * pi)

print("The Gromov-Hausdorff distance 'd' is given by the solution to the equation: cos(1 / (2d)) + cos(pi / (2d)) = 0")
print("This can be solved analytically, yielding d = (pi + 1) / (2 * pi)")
print(f"d = ({pi} + 1) / (2 * {pi})")
print(f"d = 1/2 + 1/(2*pi) = {val_1_div_2} + {val_1_div_2pi}")
print(f"The final distance is: {d}")