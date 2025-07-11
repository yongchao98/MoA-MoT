import math

# Step 1: Determine n1 and n2 from the problem's conditions.
# The condition is that u_r(n) = n/2 - 1 is the order of a differential
# equation, so it must be a positive integer.
# For n/2 - 1 to be an integer, n must be an even number.
# For n/2 - 1 > 0, we must have n/2 > 1, which means n > 2.
# The 1st and 2nd smallest positive integers n that are even and greater than 2
# are 4 and 6.
n1 = 4
n2 = 6
print(f"Step 1: Deduced from the conditions, the integers are n1 = {n1} and n2 = {n2}.")

# Step 2: Simplify the Hamiltonian H(p, q).
# The Hamiltonian is H(p, q) = 1/2 * (p^2 + q^2 - 2/n1 * sqrt((n2 - n1) / (n1 / 2)) * q^(n1 / 2)).
# We calculate the coefficient and the exponent of the q term in the potential.
exponent = n1 / 2
sqrt_term_num = n2 - n1
sqrt_term_den = n1 / 2
coeff_multiplier = math.sqrt(sqrt_term_num / sqrt_term_den)
coeff = (2 / n1) * coeff_multiplier

# The potential is V(q) = 1/2 * (q^2 - coeff * q^exponent).
# With n1=4, n2=6:
# exponent = 4 / 2 = 2
# coeff = (2/4) * sqrt((6-4)/(4/2)) = 0.5 * sqrt(2/2) = 0.5 * 1 = 0.5
# V(q) = 1/2 * (q^2 - 0.5 * q^2) = 1/2 * (0.5 * q^2) = 1/4 * q^2.
# Thus, the Hamiltonian simplifies to H = 1/2 * p^2 + 1/4 * q^2.
print("Step 2: The Hamiltonian simplifies to H = 1/2 * p^2 + 1/4 * q^2.")

# Step 3: Calculate the period of the simplified system.
# This is a harmonic oscillator H = p^2/(2m) + 1/2 * k * q^2.
# By comparing terms, we find mass m = 1 and spring constant k = 1/2.
m = 1.0
k = 0.5
# The angular frequency omega is sqrt(k/m).
omega = math.sqrt(k / m)
# The period T is 2 * pi / omega.
# T = 2 * pi / sqrt(0.5) = 2 * pi / (1/sqrt(2)) = 2 * pi * sqrt(2).
print(f"Step 3: The system is a harmonic oscillator with period T = 2 * pi * sqrt(2).")

# Step 4: Output the components of the final equation and the result.
# The final equation is T = 2 * pi * sqrt(2).
print("\nFinal Equation Details:")
val_2 = 2
val_pi = math.pi
val_sqrt2 = math.sqrt(2)
print(f"The number 2: {val_2}")
print(f"The number pi: {val_pi}")
print(f"The number sqrt(2): {val_sqrt2}")

# The final answer is T((n1-1)/n2)
final_answer = val_2 * val_pi * val_sqrt2
alpha_val = (n1 - 1) / n2
print(f"\nThe value of T({alpha_val}) is {final_answer}")