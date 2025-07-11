import sympy as sp

# 1. Define symbols and parametric equations
t = sp.Symbol('t')
x = sp.cos(t)**3
y = sp.sin(t)**3

# 2. Calculate the derivatives
dx_dt = sp.diff(x, t)
dy_dt = sp.diff(y, t)

# 3. Calculate the integrand for the arc length
# The simplify function helps to apply trigonometric identities
integrand = sp.sqrt(sp.simplify(dx_dt**2 + dy_dt**2))
print(f"The arc length integrand is: {integrand}")

# 4. Calculate the length of the arc in the first quadrant (t from 0 to pi/2)
# In this interval, sin(t) and cos(t) are non-negative, so Abs can be removed.
integrand_q1 = 3 * sp.sin(t) * sp.cos(t)
length_q1 = sp.integrate(integrand_q1, (t, 0, sp.pi/2))
print(f"The length of the arc in the first quadrant is: {length_q1}")
print(f"This matches the given arc length of 3/2.")

# 5. Determine the range of x for this arc
x_at_t0 = x.subs(t, 0)
x_at_t_pi_half = x.subs(t, sp.pi/2)
print(f"The x-coordinate at t=0 is: {x_at_t0}")
print(f"The x-coordinate at t=pi/2 is: {x_at_t_pi_half}")
print("The range of x-values for this arc is [0, 1].")

# 6. Conclude the value of 'a'
# The problem states the range of x is [0, a].
# By comparing the calculated range [0, 1] with the given range [0, a], we find a.
a = 1
print("\nEquating the calculated x-range [0, 1] with the given range [0, a] gives the final equation.")
print(f"The final equation is: a = {a}")
print(f"The value of a is: {a}")
