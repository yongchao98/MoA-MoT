import math

# Step 1 & 2: Determine Φ(z) and the integrand g(x)
# From the analysis of the fractional derivative constraints, the heat flux is:
# Φ(z) = (3/4)z - (1/π)sqrt(z)
# From the analysis of the problem structure and the limits z=1 and z=4,
# the integrand's argument is g(x) = T_1(sqrt(2)x) + T_2(x) = 3x + 1.
# This leads to the integral:
# E_total = integral from 0 to 1 of Φ(3x + 1) dx

# Step 3: Analytically evaluate the integral
# E_total = integral_0^1 [ (3/4)(3x+1) - (1/pi)*sqrt(3x+1) ] dx
# Using u-substitution with u = 3x+1, du = 3dx.
# The limits change from [0, 1] to [1, 4].
# E_total = (1/3) * integral_1^4 [ (3/4)u - (1/pi)*sqrt(u) ] du
# E_total = (1/3) * [ (3/8)u^2 - (2/(3*pi))*u^(3/2) ] from 1 to 4
# E_total = (1/3) * { [ (3/8)*16 - (2/(3*pi))*8 ] - [ (3/8)*1 - (2/(3*pi))*1 ] }
# E_total = (1/3) * { [ 6 - 16/(3*pi) ] - [ 3/8 - 2/(3*pi) ] }
# E_total = (1/3) * { (6 - 3/8) - (16/(3*pi) - 2/(3*pi)) }
# E_total = (1/3) * { 45/8 - 14/(3*pi) }
# E_total = 15/8 - 14/(9*pi)

# Step 4: Calculate the numerical value and print the components
a = 15
b = 8
c = 14
d = 9

energy_value = a / b - c / (d * math.pi)

print("The final expression for the total heat energy is E = a/b - c/(d*pi)")
print(f"The integer components of the final equation are:")
print(f"a = {a}")
print(f"b = {b}")
print(f"c = {c}")
print(f"d = {d}")
print("\nThe final equation is E = 15/8 - 14/(9*pi)")
print(f"The minimum value of the total heat energy is: {energy_value}")