import numpy as np

# The problem reduces to finding the value of Y, which is the real root of the cubic equation Y^3 - 18Y - 36 = 0.
# We define the coefficients of the polynomial P(Y) = Y^3 + 0*Y^2 - 18Y - 36.
coefficients = [1, 0, -18, -36]

# Find the roots of the polynomial equation.
roots = np.roots(coefficients)

# The result contains one real root and two complex conjugate roots. We need the real one.
Y = roots[np.isreal(roots)].real[0]

# The initial position is x(0) = 3 + Y
x0 = 3 + Y

# The time is t = 2*sqrt(3)
t_val = 2 * np.sqrt(3)

# The displacement is -t^2 / 4
displacement = -t_val**2 / 4

# The final position is x(t) = x(0) - t^2/4
xt = x0 + displacement

# We print the components of the final calculation
print(f"The equation for the particle's trajectory is x(t) = x(0) - t^2/4.")
print(f"The initial position is x(0) = 3 + Y, where Y is the real root of Y^3 - 18Y - 36 = 0.")
print(f"The numerical value of Y is: {Y}")
print(f"So, the initial position x(0) is approximately: {x0}")
print(f"At time t = 2*sqrt(3), the displacement is -t^2/4 = {displacement}")
print(f"The final position is x(2*sqrt(3)) = {x0} + ({displacement})")
print("Final Answer:")
print(f"{xt}")