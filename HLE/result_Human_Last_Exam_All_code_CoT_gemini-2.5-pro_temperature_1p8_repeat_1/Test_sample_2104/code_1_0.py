import math

# Step 1: Define n1 and n2 based on the problem's conditions.
# The conditions imply n must be a positive even integer.
# The 1st and 2nd smallest positive even integers are 2 and 4.
n1 = 2
n2 = 4

print(f"Found n1 = {n1} and n2 = {n2}")

# Step 2: Define the Hamiltonian H(p,q) by substituting n1 and n2.
# H(p,q) = 1/2 * (p^2 + q^2 - C * q^E)
# where C = (2/n1) * sqrt((n2-n1)/(n1/2)) and E = n1/2.

exp = n1 / 2
coeff_part1 = 2 / n1
coeff_part2 = math.sqrt((n2 - n1) / (n1 / 2))
coeff = coeff_part1 * coeff_part2

print("\nCalculating coefficients for the Hamiltonian:")
print(f"The exponent on q is n1/2 = {n1}/2 = {exp}")
print("The coefficient on the q term is (2/n1) * sqrt((n2-n1)/(n1/2))")
print(f"= (2/{n1}) * sqrt(({n2}-{n1})/({n1}/2))")
print(f"= {coeff_part1} * sqrt({n2 - n1}/{n1 / 2})")
print(f"= {coeff_part1} * sqrt({(n2-n1)/(n1/2)})")
print(f"= {coeff_part1} * {coeff_part2}")
print(f"= {coeff}")

print("\nThe Hamiltonian's potential energy V(q) is:")
print(f"V(q) = 1/2 * (q^2 - {coeff} * q^{exp})")
print("This describes a simple harmonic oscillator.")

# Step 3: Determine the period of this oscillator.
# The potential is V(q) = 1/2 * q^2 - (sqrt(2)/2)*q.
# This corresponds to a harmonic oscillator with mass m=1 and spring constant k=1.
# The period T_period = 2 * pi * sqrt(m/k).
period = 2 * math.pi
print(f"\nThe period of this oscillator is 2*pi, which is approximately {period}.")
print("The associated hypergeometric period function T(alpha) is therefore the constant function T(alpha) = 2*pi.")

# Step 4: Calculate the final answer.
# We need to evaluate T(alpha) at alpha = (n1-1)/n2.
alpha = (n1 - 1) / n2
print(f"\nWe need to find T(alpha) for alpha = (n1-1)/n2 = ({n1}-1)/{n2} = {alpha}.")
final_answer = period
print(f"T({alpha}) = {final_answer}")
