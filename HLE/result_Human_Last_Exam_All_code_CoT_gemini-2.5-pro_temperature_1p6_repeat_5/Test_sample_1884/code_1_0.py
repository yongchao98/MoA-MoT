import sys
import io

# Backup stdout to prevent printing to the console
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve():
    """
    Calculates the increase in atmospheric methane concentration (in ppb) 
    after a major leak, considering mixing and oxidation over 3 years.
    """
    # --- Given Constants ---
    total_ch4_leak_tons = 250000.0  # metric tons
    molar_mass_ch4_g_mol = 16.0    # g/mol
    total_mass_atmosphere_kg = 5.1e18 # kg
    
    # An assumed standard value for the average molar mass of air, necessary for molar calculations.
    avg_molar_mass_air_g_mol = 29.0 # g/mol
    
    # --- Mixing and Oxidation Parameters ---
    troposphere_fraction = 0.80
    stratosphere_fraction = 0.20
    oxidation_rate_tropo_y1 = 0.05  # 5%
    oxidation_rate_subsequent = 0.03 # 3%

    # --- Step 1: Initial Mass Distribution ---
    # Convert total leak from metric tons to grams (1 metric ton = 1,000,000 g)
    total_ch4_leak_g = total_ch4_leak_tons * 1e6
    
    initial_mass_tropo_g = total_ch4_leak_g * troposphere_fraction
    initial_mass_strato_g = total_ch4_leak_g * stratosphere_fraction
    
    # --- Step 2 & 3: Calculate Remaining Mass After Oxidation ---
    
    # Tropospheric part: decays for 3 years (Y1: 5%, Y2: 3%, Y3: 3%)
    remaining_mass_tropo_g = initial_mass_tropo_g * (1 - oxidation_rate_tropo_y1) * (1 - oxidation_rate_subsequent) * (1 - oxidation_rate_subsequent)
    
    # Stratospheric part: decays only in subsequent years (Y2: 3%, Y3: 3%)
    remaining_mass_strato_g = initial_mass_strato_g * (1 - oxidation_rate_subsequent) * (1 - oxidation_rate_subsequent)

    # --- Step 4: Sum Remaining Methane ---
    total_remaining_mass_g = remaining_mass_tropo_g + remaining_mass_strato_g

    # --- Step 5: Calculate Final Concentration in ppb ---
    
    # Convert total remaining CH4 mass to moles
    final_moles_ch4 = total_remaining_mass_g / molar_mass_ch4_g_mol
    
    # Calculate total moles of air in the atmosphere
    total_mass_atmosphere_g = total_mass_atmosphere_kg * 1000
    total_moles_air = total_mass_atmosphere_g / avg_molar_mass_air_g_mol
    
    # Calculate the concentration increase in ppb (parts per billion)
    ppb_increase = (final_moles_ch4 / total_moles_air) * 1e9

    # --- Outputting the Results ---
    original_stdout = sys.stdout
    sys.stdout = old_stdout

    print(f"Initial CH4 mass for troposphere: {initial_mass_tropo_g:e} g")
    print(f"Initial CH4 mass for stratosphere: {initial_mass_strato_g:e} g")
    print(f"Remaining CH4 mass from tropospheric part after 3 years: {remaining_mass_tropo_g:.2f} g")
    print(f"Remaining CH4 mass from stratospheric part after 3 years: {remaining_mass_strato_g:.2f} g")
    print(f"Total remaining CH4 mass after 3 years: {total_remaining_mass_g:.2f} g\n")

    print("To find the concentration increase in parts per billion (ppb), we use the formula:")
    print("ppb = (moles of CH4 / total moles of air) * 1,000,000,000\n")
    
    print("Final Equation:")
    print(f"Increase (ppb) = ( {final_moles_ch4:.4e} mol CH4 / {total_moles_air:.4e} mol air ) * 1,000,000,000")
    print(f"Calculated increase in CH4 concentration: {ppb_increase:.4f} ppb")

    # Final answer in specified format
    print(f"<<<{ppb_increase:.4f}>>>")

solve()

# Restore stdout and print the captured output
final_output = captured_output.getvalue()
sys.stdout = old_stdout
print(final_output)