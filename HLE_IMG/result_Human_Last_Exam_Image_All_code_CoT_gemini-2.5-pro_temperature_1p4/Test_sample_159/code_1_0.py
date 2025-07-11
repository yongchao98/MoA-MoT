import math

# 1. Define constants and assumptions from the problem
V_fundamental_peak = 1.0  # V
f_fundamental = 915e6      # Hz
f0 = 915e6                 # Hz (base frequency for parasitic resistance)
R0 = 50.0                  # Ohms (base parasitic resistance)
R_L = 8000.0               # Ohms (Load resistance)

# Assumption for rectifier's effective input resistance
R_rectifier = R_L / 2.0

# Harmonics and voltage drop factor
harmonics = [1, 3, 5, 7]
voltage_drop_factor = 0.9

# 2. Initialize total power variables
total_power_rectifier = 0.0
total_power_loss = 0.0

# 3. Loop through harmonics to calculate and sum powers
for n in harmonics:
    # Calculate peak voltage for the current harmonic
    V_peak = V_fundamental_peak * (voltage_drop_factor ** ((n - 1) / 2))
    
    # Calculate frequency-dependent parasitic resistance
    R_parasitic = R0 * (n**2)
    
    # Calculate total resistance for the series model
    R_total = R_parasitic + R_rectifier
    
    # Calculate power delivered to rectifier (useful) for this harmonic
    # P = V_peak^2 * R / (2 * R_total^2)
    power_rectifier_n = (V_peak**2 * R_rectifier) / (2 * R_total**2)
    
    # Calculate power lost in parasitic resistor for this harmonic
    power_loss_n = (V_peak**2 * R_parasitic) / (2 * R_total**2)
    
    # Add to totals
    total_power_rectifier += power_rectifier_n
    total_power_loss += power_loss_n

# 4. Calculate total input power and efficiency
total_input_power = total_power_rectifier + total_power_loss
efficiency = (total_power_rectifier / total_input_power)

# 5. Print the results showing the steps of the final calculation
print("--- Final Calculation ---")
print(f"Total useful power delivered to rectifier (P_useful): {total_power_rectifier:.4g} W")
print(f"Total power lost in parasitics (P_loss):         {total_power_loss:.4g} W")
print("-" * 35)
print(f"Total input power (P_in = P_useful + P_loss):   {total_input_power:.4g} W\n")

print("Efficiency (η) Calculation:")
print(f"η = P_useful / P_in")
print(f"η = {total_power_rectifier:.4g} W / {total_input_power:.4g} W")
print(f"η = {efficiency:.3f}")
print(f"η = {efficiency * 100:.1f} %\n")
print(f"The overall system efficiency is {efficiency * 100:.1f}%.")
<<<87.0>>>