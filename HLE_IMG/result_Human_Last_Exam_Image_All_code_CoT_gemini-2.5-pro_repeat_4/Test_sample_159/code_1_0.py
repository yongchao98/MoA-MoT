import math

# 1. Define all the given and assumed constants.
V_RF = 1.0  # Peak voltage of the fundamental frequency in Volts
R0 = 50.0  # Base parasitic resistance in Ohms
f0 = 915e6  # Fundamental frequency in Hz
# Assumption: The input impedance of the rectifier circuit is 50 Ohms.
R_circuit = 50.0

# 2. Define the harmonics and their peak voltages.
harmonics = [1, 3, 5, 7]
V_peaks = {}
for i, n in enumerate(harmonics):
    # V_peak for harmonic n is V_RF * (0.9)^i
    V_peaks[n] = V_RF * (0.9**i)

# 3. Loop through each harmonic to calculate its contribution to the input power.
P_in_total = 0.0
P_in_components = {}

print("Calculating Input Power Components:")
for n in harmonics:
    # Calculate frequency-dependent parasitic resistance
    R_p = R0 * (n**2)
    
    # Calculate the total impedance seen by the source
    Z_total = R_p + R_circuit
    
    # Calculate the power drawn from the source for this harmonic
    V_n = V_peaks[n]
    P_in_n = (V_n**2) / (2 * Z_total)
    P_in_components[n] = P_in_n
    P_in_total += P_in_n
    
    print(f"  For harmonic n={n}:")
    print(f"    V_peak({n}) = {V_n:.4f} V")
    print(f"    R_parasitic({n}) = {R_p:.1f} Ohms")
    print(f"    P_in({n}) = ({V_n:.4f}^2) / (2 * ({R_p:.1f} + {R_circuit:.1f})) = {P_in_n*1000:.4f} mW")

print(f"\nTotal Input Power (P_in_total) = {P_in_total*1000:.4f} mW\n")

# 4. Calculate the useful output power.
# This is the power from the fundamental (n=1) delivered to the circuit.
n1 = 1
R_p1 = R0 * (n1**2)
Z_total1 = R_p1 + R_circuit
P_out = P_in_components[n1] * R_circuit / Z_total1

print("Calculating Useful Output Power (P_out):")
print(f"  P_out = P_in(1) * R_circuit / (R_parasitic(1) + R_circuit)")
print(f"  P_out = {P_in_components[n1]*1000:.4f} mW * {R_circuit:.1f} / ({R_p1:.1f} + {R_circuit:.1f})")
print(f"  P_out = {P_out*1000:.4f} mW\n")

# 5. Calculate and print the final efficiency.
efficiency = P_out / P_in_total

print("Final Efficiency Calculation:")
print(f"η = P_out / P_in_total")
p_in_values = [p * 1000 for p in P_in_components.values()]
p_in_str = " + ".join([f"{val:.4f}" for val in p_in_values])
print(f"η = {P_out*1000:.4f} / ({p_in_str})")
print(f"η = {P_out*1000:.4f} / {P_in_total*1000:.4f}")
print(f"η = {efficiency:.4f}")
print(f"\nOverall System Efficiency = {efficiency * 100:.1f}%")
