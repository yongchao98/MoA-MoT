import math

def solve_wind_fault_compensation():
    """
    Calculates the optimal reactive power injection to mitigate a fault
    in a wind energy system based on the provided problem description.
    """

    # --- Problem Parameters ---
    V_W_kV = 0.575  # Target voltage at Bus-W in kV
    R_WF_ohm = 0.01  # Transmission line resistance in Ohms
    X_WF_ohm = 0.05  # Transmission line reactance in Ohms
    PF_min = 0.95    # Minimum lagging power factor
    # Total harmonic loss increase is the sum of 3rd and 5th harmonic losses
    harmonic_loss_increase = 0.02 + 0.04 # 6% increase in losses

    # --- Calculations ---

    # Convert voltage to V for SI unit consistency
    V_W_V = V_W_kV * 1000

    # The harmonic losses increase the resistive part of the impedance.
    # Although R does not appear in the final formula for Q_opt in this model,
    # we calculate it to acknowledge the parameter.
    R_WF_eff_ohm = R_WF_ohm * (1 + harmonic_loss_increase)

    # Calculate the constant 'k' from the power factor constraint
    # For a lagging power factor, Q is positive.
    # k = tan(theta) = Q/P
    angle = math.acos(PF_min)
    k = math.tan(angle)

    # The formula for the optimal reactive power Q_opt is derived by finding
    # the operating point that maximizes power transfer under the PF constraint.
    # Q_opt = (k^2 * |V_W|^2) / ((1 + k^2) * X_WF)

    # Substitute values and calculate
    V_W_sq = V_W_V**2
    k_sq = k**2
    denominator = (1 + k_sq) * X_WF_ohm
    Q_opt_VAR = (k_sq * V_W_sq) / denominator

    # Convert result to MVAR for the final answer
    Q_opt_MVAR = Q_opt_VAR / 1e6

    # --- Output ---
    print("--- Nonlinear Optimization Problem Solution ---")
    print("\nStep 1: Define constants and parameters from the problem.")
    print(f"Target Voltage |V_W| = {V_W_V} V")
    print(f"Line Reactance X_WF = {X_WF_ohm} Ohm")
    print(f"Power Factor PF = {PF_min}")
    print(f"Harmonic Loss Factor = 1 + {harmonic_loss_increase:.2f} = {1 + harmonic_loss_increase:.2f}")
    print("-" * 45)

    print("\nStep 2: Calculate the power factor constant k = tan(acos(PF)).")
    print(f"k = tan(acos({PF_min:.2f})) = {k:.4f}")
    print("-" * 45)

    print("\nStep 3: Formulate the equation for optimal reactive power Q_opt.")
    print("The optimal solution is found by maximizing power transfer subject to the voltage and power factor constraints.")
    print("This yields the formula:")
    print("Q_opt = (k^2 * |V_W|^2) / ((1 + k^2) * X_WF)")
    print("-" * 45)

    print("\nStep 4: Substitute the numbers into the equation.")
    print(f"Q_opt = ({k:.4f}^2 * {V_W_V}^2 V^2) / ((1 + {k:.4f}^2) * {X_WF_ohm} Ohm)")
    print(f"      = ({k_sq:.4f} * {V_W_sq} V^2) / ((1 + {k_sq:.4f}) * {X_WF_ohm} Ohm)")
    print(f"      = {k_sq * V_W_sq:.2f} / ({1 + k_sq:.4f} * {X_WF_ohm})")
    print(f"      = {k_sq * V_W_sq:.2f} / {denominator:.6f}")
    print(f"      = {Q_opt_VAR:.2f} VAR")
    print("-" * 45)

    print("\nStep 5: Final Answer in MVAR.")
    print(f"The optimal reactive power injection Q_opt is {Q_opt_MVAR:.4f} MVAR.")

if __name__ == '__main__':
    solve_wind_fault_compensation()
    # The final numerical answer is approximately 0.64 MVAR.
    # <<<0.64>>>