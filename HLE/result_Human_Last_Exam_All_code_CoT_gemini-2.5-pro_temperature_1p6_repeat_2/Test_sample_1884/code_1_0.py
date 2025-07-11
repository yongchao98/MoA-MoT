import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = buffer = io.StringIO()

def solve_methane_concentration():
    """
    Calculates the increase in atmospheric methane concentration after a major leak,
    considering atmospheric mixing and oxidation over 3 years.
    """
    # --- Step 1: Define Constants and Initial Values ---
    total_ch4_leak_tons = 250000  # metric tons
    molar_mass_ch4 = 16.0  # g/mol
    molar_mass_air = 29.0  # g/mol (average)
    mass_troposphere_kg = 4.5e18  # kg
    mass_stratosphere_kg = 0.6e18 # kg

    # --- Step 2: Calculate Moles of CH4 and Air ---
    # Convert CH4 leak to grams, then moles
    total_ch4_leak_g = total_ch4_leak_tons * 1e6
    total_moles_ch4 = total_ch4_leak_g / molar_mass_ch4

    # Distribute CH4 moles into troposphere and stratosphere
    moles_ch4_tropo_initial = total_moles_ch4 * 0.80
    moles_ch4_strato_reservoir = total_moles_ch4 * 0.20
    moles_ch4_strato_yearly_influx = moles_ch4_strato_reservoir / 3.0

    # Calculate moles of air in each layer
    moles_air_troposphere = (mass_troposphere_kg * 1000) / molar_mass_air
    moles_air_stratosphere = (mass_stratosphere_kg * 1000) / molar_mass_air
    total_moles_air = moles_air_troposphere + moles_air_stratosphere

    # --- Step 3: Track Concentration Changes (in ppb) Over 3 Years ---
    # ppb = (moles_gas / moles_air) * 1e9

    # Initial state (Year 0)
    tropospheric_increase_ppb = (moles_ch4_tropo_initial / moles_air_troposphere) * 1e9
    stratospheric_increase_ppb = 0.0

    # Stratospheric yearly influx in ppb
    stratospheric_influx_ppb = (moles_ch4_strato_yearly_influx / moles_air_stratosphere) * 1e9

    # Year 1
    # Troposphere: 5% reduction
    tropospheric_increase_ppb *= (1 - 0.05)
    # Stratosphere: First year's influx is added
    stratospheric_increase_ppb += stratospheric_influx_ppb

    # Year 2
    # Troposphere: 3% reduction
    tropospheric_increase_ppb *= (1 - 0.03)
    # Stratosphere: Second year's influx is added, then 3% reduction on the total
    stratospheric_increase_ppb += stratospheric_influx_ppb
    stratospheric_increase_ppb *= (1 - 0.03)

    # Year 3
    # Troposphere: 3% reduction
    tropospheric_increase_ppb *= (1 - 0.03)
    # Stratosphere: Third year's influx is added, then 3% reduction on the total
    stratospheric_increase_ppb += stratospheric_influx_ppb
    stratospheric_increase_ppb *= (1 - 0.03)
    
    final_tropospheric_increase_ppb = tropospheric_increase_ppb
    final_stratospheric_increase_ppb = stratospheric_increase_ppb

    # --- Step 4: Calculate Final Weighted Average Concentration Increase ---
    weight_troposphere = moles_air_troposphere / total_moles_air
    weight_stratosphere = moles_air_stratosphere / total_moles_air

    final_average_increase_ppb = (final_tropospheric_increase_ppb * weight_troposphere) + \
                                 (final_stratospheric_increase_ppb * weight_stratosphere)

    # --- Step 5: Print the Results ---
    print("Calculation of Final Weighted Average Methane Concentration Increase:")
    print("-" * 65)
    print(f"Final Tropospheric CH4 Increase: {final_tropospheric_increase_ppb:.3f} ppb")
    print(f"Weight of Troposphere (by mole fraction): {weight_troposphere:.4f}")
    print(f"Final Stratospheric CH4 Increase: {final_stratospheric_increase_ppb:.3f} ppb")
    print(f"Weight of Stratosphere (by mole fraction): {weight_stratosphere:.4f}")
    print("-" * 65)
    
    # Print the equation as requested
    term1_val = final_tropospheric_increase_ppb * weight_troposphere
    term2_val = final_stratospheric_increase_ppb * weight_stratosphere
    
    print("Final Equation:")
    print(f"({final_tropospheric_increase_ppb:.3f} ppb * {weight_troposphere:.4f}) + ({final_stratospheric_increase_ppb:.3f} ppb * {weight_stratosphere:.4f}) = {final_average_increase_ppb:.3f} ppb")
    print(f"{term1_val:.3f} + {term2_val:.3f} = {final_average_increase_ppb:.3f} ppb")
    print("-" * 65)

    print(f"Total Increase in Atmospheric Methane Concentration After 3 Years: {final_average_increase_ppb:.3f} ppb")

    # Store final numerical answer for extraction
    global final_answer
    final_answer = final_average_increase_ppb

# Run the calculation
solve_methane_concentration()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = buffer.getvalue()

# Print the captured output
print(output)

# Print the final answer in the required format
print(f"<<<{final_answer:.3f}>>>")