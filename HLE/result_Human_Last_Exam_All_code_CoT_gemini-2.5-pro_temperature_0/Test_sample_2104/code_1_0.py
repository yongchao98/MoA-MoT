import math

# Step 1: Define n1 and n2 based on the problem's conditions.
# The conditions identify n as even integers where the order of the Picard-Fuchs equation
# is reduced from the standard n-2 to n/2 - 1. This is a known property for n=4 and n=6.
n1 = 4
n2 = 6

# Step 2: Calculate the parameters for the Hamiltonian's potential V(q).
# The Hamiltonian is H(p, q) = 1/2 * (p^2 + q^2 - C * q^k)
# where C = 2/n1 * sqrt((n2 - n1) / (n1 / 2)) and k = n1 / 2.

# Calculate the coefficient C
c_numerator = n2 - n1
c_denominator = n1 / 2
c_sqrt_term = math.sqrt(c_numerator / c_denominator)
C = (2 / n1) * c_sqrt_term

# Calculate the exponent k
k = n1 / 2

# Step 3: Simplify the potential V(q) and the Hamiltonian H(p, q).
# The potential is V(q) = 1/2 * (q^2 - C * q^k).
# Since k=2, the potential simplifies to a purely quadratic form: V(q) = c * q^2.
# The Hamiltonian becomes H(p, q) = p^2/2 + c * q^2.
c = 0.5 * (1 - C)

# Step 4: Calculate the period T of the resulting system.
# The system is a simple harmonic oscillator H = p^2/(2m) + (1/2)k_s*q^2.
# By comparison, mass m=1 and spring constant k_s = 2*c.
# The angular frequency is omega = sqrt(k_s / m) = sqrt(2*c).
# The period is T = 2*pi / omega.
omega = math.sqrt(2 * c)
period = 2 * math.pi / omega

# Step 5: Evaluate T(alpha) for the given alpha.
# The problem asks for T((n1 - 1) / n2).
# This is interpreted as evaluating the period function T(E) at energy E = alpha.
# Since the period of this specific harmonic oscillator is constant, T(E) = T for any E.
alpha_numerator = n1 - 1
alpha_denominator = n2
alpha = alpha_numerator / alpha_denominator

# Step 6: Output the numbers in the final equation.
# The final equation is T((n1-1)/n2) = result.
# We will print the components that make up this equation.
print("--- Calculation Steps ---")
print(f"1. Integers from problem conditions: n1 = {n1}, n2 = {n2}")
print(f"2. Argument for T: alpha = ({n1} - 1) / {n2} = {alpha}")
print(f"3. Hamiltonian coefficient C = (2/{n1})*sqrt(({n2}-{n1})/({n1}/2)) = {C}")
print(f"4. Hamiltonian exponent k = {n1}/2 = {k}")
print(f"5. Simplified potential V(q) = (0.5 * (1 - C)) * q^2 = {c} * q^2")
print(f"6. Angular frequency omega = sqrt(2 * c) = sqrt(2 * {c}) = {omega}")
print(f"7. Period T = 2*pi/omega = {period}")

print("\n--- Final Equation ---")
# The final equation shows the substitution of all calculated values.
print(f"T(({n1} - 1) / {n2}) = 2 * pi / sqrt(2 * (0.5 * (1 - (2/{n1}) * sqrt(({n2}-{n1})/({n1}/2)))))")
print(f"T({alpha}) = {period}")