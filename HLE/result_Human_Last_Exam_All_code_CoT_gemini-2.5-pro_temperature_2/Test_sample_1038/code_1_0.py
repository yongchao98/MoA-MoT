import math

def calculate_current_ratio(v_gs, vt_normal, vt_biased, n_factor, v_thermal_mV):
    """
    Calculates the ratio of subthreshold current for a body-biased transistor
    compared to a normal transistor.

    The subthreshold current is proportional to exp((Vgs - Vt) / (n * V_thermal)).
    The ratio of (current_biased / current_normal) shows the change in leakage.
    The resistance is inversely proportional to this ratio.
    """
    
    # Convert thermal voltage from mV to V
    v_thermal = v_thermal_mV / 1000.0

    # Calculate the exponent term for the normal case
    exponent_normal = (v_gs - vt_normal) / (n_factor * v_thermal)
    
    # Calculate the exponent term for the body-biased case (increased Vt)
    exponent_biased = (v_gs - vt_biased) / (n_factor * v_thermal)

    # The currents are proportional to math.exp() of these values.
    # The ratio of currents is exp(exponent_biased) / exp(exponent_normal)
    current_ratio = math.exp(exponent_biased - exponent_normal)
    
    # Resistance is inversely proportional to current
    resistance_increase_factor = 1 / current_ratio
    
    return current_ratio, resistance_increase_factor

# --- Circuit Parameters ---
# Nominal Threshold Voltage (Volts)
vt_normal = 0.45

# We need a Vgs that places the transistor in subthreshold, e.g., slightly below Vt.
# The exact value isn't critical, as we are looking at the ratio.
# Let's assume the gate is biased 100 mV below the nominal threshold.
v_gs = vt_normal - 0.100

# Subthreshold slope factor (dimensionless)
n_factor = 1.5

# Thermal voltage at room temperature (mV)
v_thermal_mV = 26 

# --- Option C Strategy: Body Biasing ---
# Body biasing can effectively increase the threshold voltage. Let's assume
# it increases Vt by 150 mV during the 'operate' phase.
vt_increase = 0.150
vt_biased = vt_normal + vt_increase

# --- Calculation ---
current_ratio, resistance_factor = calculate_current_ratio(
    v_gs=v_gs,
    vt_normal=vt_normal,
    vt_biased=vt_biased,
    n_factor=n_factor,
    v_thermal_mV=v_thermal_mV
)

# --- Output the results ---
print("Analysis of Body Biasing Strategy (Option C):")
print("-" * 50)
print(f"Assumed Gate-Source Voltage (Vgs): {v_gs:.2f} V")
print(f"Nominal Threshold Voltage (Vt): {vt_normal:.2f} V")
print(f"Body-Biased Threshold Voltage (Vt_biased): {vt_biased:.2f} V (Vt + {vt_increase:.3f} V)")
print(f"Subthreshold Slope Factor (n): {n_factor}")
print("-" * 50)
print(f"The subthreshold current with body biasing is {current_ratio:.4f} times the original current.")
print(f"This means the effective resistance is increased by a factor of approximately {resistance_factor:.1f}.")
print("\nConclusion: Body biasing offers a powerful method to drastically increase resistance")
print("during the operate phase without compromising the reset phase, making it the most effective strategy.")

<<<C>>>