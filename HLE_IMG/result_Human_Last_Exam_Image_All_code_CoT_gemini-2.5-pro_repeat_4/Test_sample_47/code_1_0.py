from fractions import Fraction

# Step 1 & 2: Define component values and known voltages from the circuit diagram.
R76 = Fraction(76)
R8 = Fraction(8)
R14 = Fraction(14)
R11 = Fraction(11)
R29 = Fraction(29)
V_source = Fraction(41)

# Let's define the nodes based on a grid interpretation.
# V_ml (mid-left node) is connected to the positive terminal of the battery.
V_ml = V_source
# V_tr (top-right node) and V_br (bottom-right node) are connected to ground.
V_tr = Fraction(0)
V_br = Fraction(0)

# Step 3 & 4: Solve for unknown node voltages V_tl (top-left) and V_mr (mid-right).

# KCL at node V_tl: (V_tl - V_tr)/R76 + (V_tl - V_ml)/R8 = 0
# V_tl * (1/R76 + 1/R8) = V_ml/R8
V_tl = (V_ml / R8) / (1/R76 + 1/R8)

# KCL at node V_mr: (V_mr - V_ml)/R14 + (V_mr - V_tr)/R11 + (V_mr - V_br)/R29 = 0
# V_mr * (1/R14 + 1/R11 + 1/R29) = V_ml/R14
V_mr = (V_ml / R14) / (1/R14 + 1/R11 + 1/R29)

# Step 5: Calculate the current I.
# I is the total current leaving the positive node V_ml.
# It splits into two paths: through R8 and through R14.
I_R8 = (V_ml - V_tl) / R8
I_R14 = (V_ml - V_mr) / R14
I_total = I_R8 + I_R14

# Print the detailed calculation for the user
print("1. Find unknown node voltages:")
print(f"   Voltage at node V_tl = (V_source / R8) / (1/R76 + 1/R8)")
print(f"   V_tl = ({V_source} / {R8}) / (1/{R76} + 1/{R8}) = {V_tl.limit_denominator()} V")
print(f"   Voltage at node V_mr = (V_source / R14) / (1/R14 + 1/R11 + 1/R29)")
print(f"   V_mr = ({V_source} / {R14}) / (1/{R14} + 1/{R11} + 1/{R29}) = {V_mr.limit_denominator()} V")
print("\n2. Calculate the current I as the sum of currents leaving the source node:")
print(f"   I = I_R8 + I_R14")
print(f"   I = (V_source - V_tl) / R8 + (V_source - V_mr) / R14")
# Note: Numerators are calculated first to show the voltage differences
v_diff_r8 = V_source - V_tl
v_diff_r14 = V_source - V_mr
print(f"   I = ({v_diff_r8.limit_denominator()}) / {R8} + ({v_diff_r14.limit_denominator()}) / {R14}")
print(f"   I = {I_R8.limit_denominator()} + {I_R14.limit_denominator()}")
print(f"   I = {I_total.limit_denominator()} A")
print(f"\nFinal Answer (as a decimal): {float(I_total):.3f} A")