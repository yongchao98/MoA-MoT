import math

# The problem is to find the integer part of 10^4 * L, where L is the limit of F(N)/ln(N).
# L is determined by the growth rates of two sequences of solutions.

# The growth rates are the dominant roots of the characteristic equations
# r^2 - 13*r + 1 = 0 and r^2 - 5*r + 1 = 0.

# Dominant root for the first sequence
phi_A = (13 + math.sqrt(165)) / 2

# Dominant root for the second sequence
phi_B = (5 + math.sqrt(21)) / 2

# The limit L = 2/ln(phi_A) + 2/ln(phi_B)
# The factor of 2 comes from pairs (a,b) and (b,a) for each sequence.
term1 = 2 / math.log(phi_A)
term2 = 2 / math.log(phi_B)
L = term1 + term2

# We need to calculate the integer part of 10^4 * L.
result = 10000 * L

print(f"The value of phi_A is (13 + sqrt(165))/2 = {phi_A}")
print(f"The value of phi_B is (5 + sqrt(21))/2 = {phi_B}")
print(f"The limit L is 2/ln(phi_A) + 2/ln(phi_B) = {term1} + {term2} = {L}")
print(f"The final value is 10^4 * L = {result}")
print(f"The integer part of the final value is {int(result)}")