import numpy as np
from scipy.optimize import minimize

def solve_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # Step 1: Define constants and harmonic properties
    V_RF = 1.0  # V, amplitude of the fundamental
    R_L = 8000  # Ohms
    R_0 = 50    # Ohms
    f_0 = 915e6 # Hz

    harmonics = [1, 3, 5, 7]
    voltages = []
    v_n = V_RF
    for i in range(len(harmonics)):
        voltages.append(v_n)
        v_n *= 0.9

    print("--- Calculating Harmonic Properties ---")
    for i, n in enumerate(harmonics):
        print(f"Harmonic {n}: Voltage amplitude V_{n} = {voltages[i]:.3f} V")
    print("-" * 20)

    # Step 2: Calculate Output Power (P_out)
    # Define the input voltage function to find its peak.
    # We minimize the negative of the function to find the maximum.
    def v_input(x):
        signal = 0
        for i, n in enumerate(harmonics):
            signal += voltages[i] * np.sin(n * x)
        return -signal

    # Numerically find the peak voltage
    # Initial guess for omega*t, peak should be between 0 and pi/2
    res = minimize(v_input, x0=0.5, bounds=[(0, np.pi)])
    V_peak = -res.fun

    # Calculate output power
    P_out = V_peak**2 / R_L

    print("--- Calculating Output Power (P_out) ---")
    print(f"Peak input voltage V_peak = {V_peak:.4f} V")
    print(f"Assuming V_DC = V_peak, output power P_out = ({V_peak:.4f})^2 / {R_L} = {P_out:.8f} W")
    print("-" * 20)

    # Step 3: Calculate Power Loss (P_loss)
    print("--- Calculating Parasitic Power Loss (P_loss) ---")
    P_loss_total = 0
    parasitic_resistances = []
    power_losses = []

    for i, n in enumerate(harmonics):
        # Calculate parasitic resistance for the n-th harmonic
        f_n = n * f_0
        R_parasitic_n = R_0 * (f_n / f_0)**2
        parasitic_resistances.append(R_parasitic_n)

        # Calculate power loss for the n-th harmonic
        P_loss_n = voltages[i]**2 / (2 * R_parasitic_n)
        power_losses.append(P_loss_n)

        print(f"Harmonic {n}:")
        print(f"  Parasitic Resistance R_parasitic({n}f_0) = {R_parasitic_n:.1f} Ohms")
        print(f"  Power Loss P_loss,{n} = ({voltages[i]:.3f})^2 / (2 * {R_parasitic_n:.1f}) = {P_loss_n:.8f} W")
        P_loss_total += P_loss_n

    print(f"\nTotal Parasitic Power Loss P_loss_total = {P_loss_total:.8f} W")
    print("-" * 20)

    # Step 4: Calculate Total Input Power (P_in)
    P_in = P_out + P_loss_total
    print("--- Calculating Total Input Power (P_in) ---")
    print(f"P_in = P_out + P_loss_total = {P_out:.8f} + {P_loss_total:.8f} = {P_in:.8f} W")
    print("-" * 20)

    # Step 5: Calculate Efficiency (η)
    efficiency = P_out / P_in
    print("--- Calculating Overall System Efficiency (η) ---")
    print(f"Efficiency η = P_out / P_in")
    print(f"η = {P_out:.8f} / ({P_out:.8f} + {P_loss_total:.8f})")
    print(f"η = {P_out:.8f} / {P_in:.8f} = {efficiency:.4f}")
    print(f"Efficiency = {efficiency * 100:.2f}%")
    
    return efficiency

# Run the calculation
final_efficiency = solve_efficiency()
# The final answer will be printed by the function above, and we also provide it in the required format.
# The value is 4.38%
final_answer_val = final_efficiency * 100
print(f"\n<<<η = {final_answer_val:.2f}%>>>")
