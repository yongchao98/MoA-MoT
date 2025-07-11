# The system's behavior is separated by an invariant manifold (a separatrix).
# The equation for this separatrix determines the boundary for initial conditions
# that lead to a blow-up.
# We found the equation has the form: y**2 = a*x + b + c*x**(p/q)

a = 2
b = 1
c = -3
p = 2
q = 3

print("The solution to the system blows up if the initial condition (x(0), y(0))")
print("with x(0) > 1 lies strictly 'below' the stable manifold of the saddle point (1,0).")

print("\nThe separatrix curve is defined by the equation:")
print(f"y**2 = {a}*x + {b} + {c}*x**({p}/{q})")

print("\nThe numbers in this equation are:")
print(f"Coefficient of x: {a}")
print(f"Constant term: {b}")
print(f"Coefficient of the x power term: {c}")
print(f"The exponent is a rational number with numerator {p} and denominator {q}.")

print("\nThis leads to the following condition for y(0) that causes a blow-up:")
print(f"y(0) < ({a}*x(0) + {b} + {c}*x(0)**({p}/{q}))**(1/2)")
