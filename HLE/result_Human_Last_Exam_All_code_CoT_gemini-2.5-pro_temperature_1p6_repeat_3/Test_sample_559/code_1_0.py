import sympy

# The system of differential equations is:
# d'(t) = 2*d(t)**2 + (-3*u(t) + 5*u(t)**2)*d(t) - u(t)*(1-u(t))*u(t)**2
# u'(t) = (u(t)-1)*u(t)**2

# Step 1: Find equilibrium points by setting d'(t)=0 and u'(t)=0.
# u'(t) = 0 gives u=0 or u=1.
# If u=0, d'(t)=0 gives 2*d**2 = 0, so d=0. Point: (0, 0).
# If u=1, d'(t)=0 gives 2*d**2 + 2*d = 0, so d=0 or d=-1. Points: (1, 0) and (1, -1).

# Step 2: Classify the equilibrium points. The point (1, -1) is a saddle point.
# A separatrix is a trajectory that separates regions of different qualitative behavior.
# The stable and unstable manifolds of a saddle point are separatrices.

# Step 3 & 4: Find the equation for the unstable manifold of the saddle point (1,-1).
# We search for an invariant curve of the form d = f(u).
# By analyzing the system, we can find a function H(u,d) = d + u**2.
# The time derivative of H is dH/dt = (2*d + 3*u**2 - 3*u) * H.
# This means that if a trajectory starts on the curve H=0, it stays on that curve.
# The equation H(u,d) = 0 gives d + u**2 = 0, or d = -u**2.
# This curve passes through the saddle point (1, -1) and the equilibrium (0, 0).
# This curve is the unstable manifold of the saddle point and is a separatrix of the system.

# The following code prints the equation of this separatrix.
# The equation is of the form: d = c * u^p

# Define the variables and parameters for the equation of the separatrix
dependent_var = "d"
independent_var = "u"
coefficient = -1
power = 2

# Print the final equation
print("The equation of the separatrix is:")
# We use string formatting to construct and print the equation.
# The numbers in the equation, -1 and 2, are explicitly handled by variables.
equation = f"{dependent_var} = {coefficient} * {independent_var}**{power}"
print(equation)
