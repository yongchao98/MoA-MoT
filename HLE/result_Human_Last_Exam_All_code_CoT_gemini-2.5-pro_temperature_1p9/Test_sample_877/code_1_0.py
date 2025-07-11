import numpy as np

# The function h(x) is derived from the separatrix equation of the system.
# h(x) = a^2, where 'a' and 'x=b' lie on the separatrix.
# The derived expression is h(x) = 4*x**2 - 6*x + 2 + 2*x*ln(2*x).

# Let's define the coefficients of the equation for clarity.
c1 = 4
c2 = -6
c3 = 2
c4 = 2  # coefficient of the x*ln(..) term
c5 = 2  # coefficient of x inside ln(..)

# We print the equation for h(x), showing each number as requested.
# This formula defines the boundary for the initial condition a(0) in terms of b(0)=x.
print("The function h(x) is determined from the separatrix of the system's phase portrait.")
print("The equation for h(x) is:")
print(f"h(x) = {c1}*x**2 + {c2}*x + {c3} + {c4}*x*ln({c5}*x)")