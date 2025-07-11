import cmath
import math

# --- 1. Define Constants based on the problem statement ---
V_W = 575.0         # Target voltage at Bus-W in Volts
R_WF = 0.01         # Transmission line resistance in Ohms
X_WF = 0.05         # Transmission line reactance in Ohms
R_F = 0.1           # Fault resistance in Ohms
PF_min = 0.95       # Minimum lagging power factor
harmonic_loss_factor = 1.06 # 1 + 2% + 4%
Q_max_MVAR = 10.0   # Maximum reactive power capacity in MVAR

# --- 2. Formulate the optimization problem ---

# The goal is to minimize Q_W subject to:
# a) A physics-based constraint g(P_W, Q_W) = 0
# b) A power factor constraint Q_W / P_W <= tan(acos(PF_min))

# The optimal solution lies on the boundary of the PF constraint,
# where PF = PF_min.
# Let k = tan(acos(PF_min)), so Q_W = k * P_W.
k = math.tan(math.acos(PF_min))

# --- 3. Derive and solve the quadratic equation for P_W ---

# The physics constraint can be simplified to a quadratic equation for P_W:
# c2 * P_W^2 + c1 * P_W + c0 = 0
# where the coefficients c2, c1, and c0 depend on the system parameters.

V2 = V_W**2

# Coefficients for the intermediate terms
A = (R_WF + k * X_WF) / V2
B = (X_WF - k * R_WF) / V2
C = harmonic_loss_factor * (1 + k**2) / V2 * R_WF

# Coefficients of the quadratic equation a*x^2 + b*x + c = 0 (using c2, c1, c0 for clarity)
c2 = (V2 * (A**2 + B**2) / R_F) + C
c1 = -(1 + (2 * A * V2) / R_F)
c0 = V2 / R_F

# --- 4. Print the equation and solve for P_W ---

print("The problem reduces to solving the following quadratic equation for active power P (in Watts):")
print(f"({c2:.6e}) * P^2 + ({c1:.6f}) * P + ({c0:.6f}) = 0")
print("-" * 30)

# Calculate the discriminant
discriminant = c1**2 - 4 * c2 * c0

if discriminant < 0:
    print("No real solutions for power exist. The system cannot be stabilized under these conditions.")
else:
    # Solve for the two possible values of P_W
    P_sol1 = (-c1 - math.sqrt(discriminant)) / (2 * c2)
    P_sol2 = (-c1 + math.sqrt(discriminant)) / (2 * c2)

    # In power systems, the lower power solution is typically the stable operating point.
    P_opt_W = min(P_sol1, P_sol2)

    # --- 5. Calculate the optimal reactive power Q_opt ---
    # We assume the wind farm has a unity power factor, so Q_comp = Q_W.
    Q_opt_VAR = k * P_opt_W
    Q_opt_MVAR = Q_opt_VAR / 1e6

    # Convert to MVA for checking capacity
    P_opt_MW = P_opt_W / 1e6

    print(f"Solved for two possible active power values: {P_sol1/1e6:.3f} MW and {P_sol2/1e6:.3f} MW.")
    print(f"Choosing the lower, stable power solution: {P_opt_MW:.3f} MW.")
    
    # Check against the compensator's capacity
    if Q_opt_MVAR > Q_max_MVAR:
        print("\nWarning: The required reactive power exceeds the maximum capacity of the device.")
        print(f"Required Q = {Q_opt_MVAR:.3f} MVAR, Max Q = {Q_max_MVAR:.3f} MVAR")
    else:
        print("\nThe optimal reactive power injection required from the compensating device is:")
        print(f"Q_opt = {Q_opt_MVAR:.4f} MVAR")
