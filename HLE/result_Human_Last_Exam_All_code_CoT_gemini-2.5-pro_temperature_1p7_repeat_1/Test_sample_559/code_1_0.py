# The system of differential equations is:
# d'(t) = 2*d(t)**2 + (-3*u(t) + 5*u(t)**2)*d(t) - u(t)*(1 - u(t))*u(t)**2
# u'(t) = (u(t) - 1)*u(t)**2

# Through analysis of the system's dynamics, we have identified the equilibrium points
# and their stability. The point (1, -1) is a saddle point. The separatrix
# is the unstable manifold of this saddle point.

# We found that the separatrix is described by a parabolic equation d = a*u**2 + b*u + c.
# The coefficients for this parabola were found to be:
a = -1
b = 0
c = 0

# We print the equation of the separatrix, including each coefficient as requested.
print("The equation of the separatrix is a curve in the (u, d) phase plane.")
print(f"The equation is: d = ({a}) * u^2 + ({b}) * u + ({c})")
print("Simplified, the equation is: d = -u^2")
