from sympy import symbols, Eq, solve

# Define the unknown node voltages
Vs_neg, V1, V2 = symbols('Vs_neg V1 V2')

# Resistor and Voltage values from the diagram
R1 = 76
R2 = 8
R3 = 14
R4 = 11
R5 = 29
V_source = 41

# Equation 1 from KCL at node V1
# (V1 - Vs_neg)/R2 + (V1 - V2)/R3 = 0
# Simplified: -7*Vs_neg + 11*V1 - 4*V2 = 0
eq1 = Eq(-7*Vs_neg + 11*V1 - 4*V2, 0)

# Equation 2 from KCL at node V2
# (V2 - V1)/R3 + V2/R4 + (V2 - (Vs_neg + V_source))/R5 = 0
# Simplified: -154*Vs_neg - 319*V1 + 879*V2 = -6314
eq2 = Eq(-154*Vs_neg - 319*V1 + 879*V2, -6314)

# Equation 3 from KCL for the voltage source
# Vs_neg/R1 + (Vs_neg - V1)/R2 = (Vs_neg + V_source - V2)/R5
# Simplified: 457*Vs_neg - 551*V1 + 152*V2 = 6232
eq3 = Eq(457*Vs_neg - 551*V1 + 152*V2, 6232)

# Solve the system of linear equations
solution = solve((eq1, eq2, eq3), (Vs_neg, V1, V2))

# Extract the voltage values
v1_val = solution[V1]
v2_val = solution[V2]

# Calculate the current I
# I = (V1 - V2) / R3
current_I = (v1_val - v2_val) / R3

print(f"Solved Node Voltages:")
print(f"V1 = {v1_val} V")
print(f"V2 = {v2_val} V")
print(f"\nCalculation for current I:")
print(f"I = (V1 - V2) / R3")
print(f"I = ({v1_val} - {v2_val}) / {R3}")
print(f"I = {current_I} A")