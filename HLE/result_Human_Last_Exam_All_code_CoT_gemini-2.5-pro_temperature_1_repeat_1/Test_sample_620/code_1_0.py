import math

def calculate_rate_constant(A, Ea, R, T_celsius):
    """
    Calculates the reaction rate constant 'k' using the Arrhenius equation.
    k = A * exp(-Ea / (R * T))
    
    Args:
        A (float): Pre-exponential factor
        Ea (float): Activation energy in J/mol
        R (float): Universal gas constant in J/(mol*K)
        T_celsius (float): Temperature in Celsius
        
    Returns:
        float: The calculated rate constant k.
    """
    # Convert temperature from Celsius to Kelvin
    T_kelvin = T_celsius + 273.15
    
    # Calculate the rate constant k
    k = A * math.exp(-Ea / (R * T_kelvin))
    
    return k, T_kelvin

# --- Parameters for the Arrhenius Equation ---
# These are typical values for an enzyme-catalyzed reaction.
A = 1e10      # Pre-exponential factor (s^-1)
Ea = 50000    # Activation energy in Joules/mol (50 kJ/mol)
R = 8.314     # Universal gas constant in J/(mol*K)

# --- Temperatures to Compare ---
T_cold_C = 0.0   # On ice
T_warm_C = 25.0  # Room temperature

# --- Calculations ---
k_cold, T_cold_K = calculate_rate_constant(A, Ea, R, T_cold_C)
k_warm, T_warm_K = calculate_rate_constant(A, Ea, R, T_warm_C)

# --- Output the results ---
print("The Arrhenius equation describes the effect of temperature on the reaction rate constant 'k'.")
print("k = A * exp(-Ea / (R * T))\n")

print(f"Analysis at low temperature ({T_cold_C}°C):")
print(f"k_cold = {A} * exp(-{Ea} / ({R} * {T_cold_K:.2f}))")
print(f"Calculated rate constant on ice (k_cold): {k_cold:.4f} s⁻¹\n")

print(f"Analysis at warm temperature ({T_warm_C}°C):")
print(f"k_warm = {A} * exp(-{Ea} / ({R} * {T_warm_K:.2f}))")
print(f"Calculated rate constant at room temp (k_warm): {k_warm:.4f} s⁻¹\n")

if k_cold > 0:
    rate_increase_factor = k_warm / k_cold
    print(f"Conclusion: Increasing the temperature from {T_cold_C}°C to {T_warm_C}°C increases the reaction rate by a factor of {rate_increase_factor:.1f}.")
    print("This demonstrates that increasing the temperature is a critical step to achieve a measurable and linear reaction rate.")
else:
    print("The rate at the cold temperature is effectively zero, highlighting the need to increase temperature.")
