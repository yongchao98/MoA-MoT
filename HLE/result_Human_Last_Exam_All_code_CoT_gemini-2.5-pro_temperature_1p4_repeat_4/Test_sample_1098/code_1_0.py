# This script models the initial and final states of the gas in the
# Maxwell's demon experiment, highlighting the role of temperature.

# Plan:
# 1. Define the initial conditions: two equal compartments with equal amounts of gas.
# 2. Set a non-zero temperature, which is the required parameter for molecular motion.
# 3. Use the Ideal Gas Law (PV = nRT) to calculate the initial pressure in both compartments.
# 4. Define the final condition: all gas has been trapped in one compartment.
# 5. Calculate the final pressure in the full compartment, assuming temperature remains constant.
# 6. Print the results, including the final pressure calculation step-by-step.

# --- Parameters and Constants ---
# Experimental parameters
n_initial_moles = 1.0 # Initial moles of gas in each compartment
volume_m3 = 1.0       # Volume of each compartment in cubic meters
temperature_K = 300   # Temperature in Kelvin. THIS IS THE REQUIRED PARAMETER.

# Physical constant
R_gas_constant = 8.314 # Ideal gas constant in J/(mol·K)

print("--- Maxwell's Demon Apparatus State ---")

# --- Initial State Calculation ---
# Pressure (P) = (n * R * T) / V
if temperature_K > 0:
    initial_pressure_Pa = (n_initial_moles * R_gas_constant * temperature_K) / volume_m3
    print(f"Required parameter Temperature is set to: {temperature_K} K\n")
    print("--- Initial State ---")
    print(f"Each compartment has {n_initial_moles} moles of gas at {initial_pressure_Pa:.2f} Pascals.\n")

    # --- Final State Calculation ---
    # All gas moves to one side. The number of moles in that side doubles.
    final_moles = n_initial_moles * 2
    final_pressure_Pa = (final_moles * R_gas_constant * temperature_K) / volume_m3

    print("--- Final State ---")
    print("After the 'demon' operates, all gas is trapped in one compartment.")
    print(f"The final pressure in that compartment is {final_pressure_Pa:.2f} Pascals.\n")
    
    # --- Final Equation Output ---
    print("Final pressure calculation:")
    print(f"P_final = (n_final * R * T) / V")
    print(f"P_final = ({final_moles} * {R_gas_constant} * {temperature_K}) / {volume_m3}")
else:
    print("Temperature must be above absolute zero for gas molecules to have motion.")
    print("Without motion, the 'demon' cannot operate.")
