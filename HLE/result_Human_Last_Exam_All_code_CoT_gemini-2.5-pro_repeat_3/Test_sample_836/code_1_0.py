import sympy
from sympy import Function, Integral, exp, I, oo, pprint

# Step 1: Define the symbolic variables and functions using sympy.
# We define the variables for space (r), time (t), wavevector (k),
# and angular frequency (omega). For simplicity, we consider one spatial dimension.
r, t, k, omega = sympy.symbols('r t k omega', real=True)

# We define the generalized pair correlation function as a symbolic function g(r, t).
# This function describes the correlation of detecting particles at different positions and times.
g = Function('g')(r, t)

# Step 2: Define the space-time double Fourier Transform.
# The transform converts the function from the space-time domain (r, t)
# to the wavevector-frequency domain (k, omega).
integrand = g * exp(-I * (k * r - omega * t))
fourier_transform_expr = Integral(integrand, (r, -oo, oo), (t, -oo, oo))

# Step 3: Print the symbolic expression and state the common name of the result.
print("The mathematical expression for the space-time double Fourier transform of the")
print("generalized pair correlation function g(r, t) can be represented as:")
print("-" * 70)
pprint(fourier_transform_expr, use_unicode=True)
print("-" * 70)
print("\nIn the nuclear criticality community, this resulting quantity is known as the:")
print(">>> Power Spectral Density (PSD) <<<\n")

# Step 4: Extract and print the numerical coefficients from the transform equation.
# The core of the transform is the exponential term: exp(-I * (k*r - omega*t)).
# Expanding the exponent gives: -I*k*r + I*omega*t.
# The numerical coefficients for the 'I*k*r' and 'I*omega*t' terms are -1 and 1.
exponent_term = (-I * (k * r - omega * t)).expand()
coeff_omega_t = exponent_term.coeff(I*omega*t)
coeff_k_r = exponent_term.coeff(I*k*r)

print("As requested, printing the numerical coefficients from the exponent of the transform equation [-I*(k*r - omega*t)]:")
print(f"The coefficient for the 'I*omega*t' term is: {coeff_omega_t}")
print(f"The coefficient for the 'I*k*r' term is: {coeff_k_r}")