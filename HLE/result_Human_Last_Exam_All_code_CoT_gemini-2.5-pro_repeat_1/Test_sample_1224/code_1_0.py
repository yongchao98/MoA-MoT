import math

# This script simulates the calculation performed on the Titan computer.
# It uses approximations that adhere to the 4-bit (0-15) integer constraint
# for numerators and denominators.

# --- Helper class to represent Titan's fractional numbers ---
class TitanNumber:
    """
    A class to represent a number in Titan's fractional system.
    It can be a simple fraction or an expression (sum of fractions).
    For this simulation, we'll just store the final float value
    but print the representation used.
    """
    def __init__(self, value, representation):
        self.value = float(value)
        self.representation = str(representation)

    def __str__(self):
        return self.representation

# --- 1. Define constants and parameters with 4-bit approximations ---
# Gravitational constant G = 6.674e-11 N(m/kg)^2
# Approximation: 2/3 * 10^-10
G = TitanNumber(2/3 * 10**-10, "(2/3)e-10")

# Pi
# Approximation: 3/1
PI = TitanNumber(3.0, "3/1")

# Pandora's shell density rho = 300 kg/m^3
# Representation: 3/1 * 10^2
RHO_SHELL = TitanNumber(300.0, "3/1e2")

# Pandora's radius r (simplified to a sphere) = 2000 km
# Representation: 2/1 * 10^6 m
R = TitanNumber(2.0e6, "2/1e6")

# Drop height d = 5000 m
# Representation: 5/1 * 10^3 m
D = TitanNumber(5000.0, "5/1e3")

# Constant 2
TWO = TitanNumber(2.0, "2/1")

# --- 2. Calculate gravitational acceleration 'g' ---
# Simplified model: g = (4/3) * pi * G * rho * r
# The calculation on Titan would involve multiplying the expression terms.
# The fractional part product is 4/3 * 3 * 2/3 * 3 * 2 = 16.
# On Titan, 16 is stored as an expression, e.g., (10/1 + 6/1).
# The exponent part product is 10^-10 * 10^2 * 10^6 = 10^-2.
g_value = (4/3) * PI.value * G.value * RHO_SHELL.value * R.value
# The representation of g is based on the result of the Titan calculation.
G_TITAN = TitanNumber(0.16, "(10/1 + 6/1)e-2")

# --- 3. Calculate landing time 't' ---
# Formula: t = sqrt(2 * d / g)
# First, calculate t^2 = 2 * d / g
t_squared_value = 2 * D.value / G_TITAN.value
# Representation on Titan: 62500 = 5/8 * 10^5
T_SQUARED_TITAN = TitanNumber(t_squared_value, "5/8e5")

# Now, find the square root. sqrt(62500) = 250.
# The number 250 can be represented as an expression: 200 + 50
# -> 2/1 * 10^2 + 5/1 * 10^1
t_final_value = math.sqrt(t_squared_value)
T_FINAL_TITAN = TitanNumber(t_final_value, "(2/1e2 + 5/1e1)")

# --- 4. Print the final equation with the computed values ---
# The final result is an equation showing how t is calculated from the other Titan numbers.
print("Final Equation on Titan:")
print(f"t = sqrt(2 * d / g)")
print(" ")
print("Calculated values:")
print(f"t = {T_FINAL_TITAN}")
print(f"2 = {TWO}")
print(f"d = {D}")
print(f"g = {G_TITAN}")
print(" ")
print("Final equation with numbers:")
print(f"{T_FINAL_TITAN} = sqrt({TWO} * {D} / {G_TITAN})")
