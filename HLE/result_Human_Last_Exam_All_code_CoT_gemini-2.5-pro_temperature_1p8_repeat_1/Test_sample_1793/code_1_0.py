import math

def calculate_tb155_yield():
    """
    Calculates the thick target yield of Tb-155 from proton irradiation of Gd2O3.
    """

    # --- 1. Define Constants and Input Parameters ---
    E_proton_initial_MeV = 15.0  # Initial proton energy in MeV
    E_proton_final_MeV = 12.0   # Final proton energy in MeV
    I_current_uA = 20.0        # Proton current in microamperes
    t_irradiation_hr = 4.0     # Irradiation time in hours
    T_half_days = 5.32         # Half-life of Tb-155 in days
    
    # Molar masses (g/mol)
    M_Gd = 157.25
    M_O = 15.999
    
    # Cross-sections (mb)
    sigma_mb = {
        16: 194.17,
        15: 182.82,
        14: 172.16,
        13: 163.3,
        12: 150.48,
    }

    # Physical constants
    N_A = 6.02214076e23  # Avogadro's number (atoms/mol)
    e_charge = 1.60217663e-19 # Elementary charge in Coulombs
    
    print("--- Calculation of Tb-155 Yield ---")
    print("\nStep 1: Define Initial Parameters")
    print(f"Initial Proton Energy: {E_proton_initial_MeV} MeV")
    print(f"Final Proton Energy: {E_proton_final_MeV} MeV")
    print(f"Proton Current: {I_current_uA} µA")
    print(f"Irradiation Time: {t_irradiation_hr} hours")
    print(f"Tb-155 Half-life: {T_half_days} days")
    
    # --- 2. Calculate Target and Beam Properties ---
    print("\nStep 2: Calculate Target and Beam Properties")
    
    # Molar Mass of Gd2O3
    M_Gd2O3 = 2 * M_Gd + 3 * M_O
    print(f"Molar Mass of Gd2O3 = 2 * {M_Gd} + 3 * {M_O} = {M_Gd2O3:.2f} g/mol")
    
    # Number of Gd atoms per gram of Gd2O3
    N_atoms_per_gram = (N_A * 2) / M_Gd2O3
    print(f"Number of Gd Atoms per gram = (2 * {N_A:.4e}) / {M_Gd2O3:.2f} = {N_atoms_per_gram:.4e} atoms/g")
    
    # Beam current in protons per second
    I_protons_per_sec = (I_current_uA * 1e-6) / e_charge
    print(f"Beam Current in Protons/sec = ({I_current_uA}e-6 C/s) / ({e_charge:.4e} C/proton) = {I_protons_per_sec:.4e} protons/s")
    
    # --- 3. Determine Target Thickness using Stopping Range Polynomial ---
    def calc_stopping_range(X):
        """Calculates stopping range Y (g/cm^2) for a proton of energy X (MeV)."""
        Y = (-0.00001208736486811230 * X**3 + 
             0.00194595770392697000 * X**2 + 
             0.00794283377547150000 * X - 
             0.00360695486492614000)
        return Y

    print("\nStep 3: Calculate Production Rate (R) Slice by Slice")
    total_production_rate = 0.0
    energy_levels = [15, 14, 13, 12]

    for i in range(len(energy_levels) - 1):
        E_start = float(energy_levels[i])
        E_end = float(energy_levels[i+1])
        
        print(f"\n--- Processing slice from {E_start} MeV to {E_end} MeV ---")

        # Calculate thickness of this slice
        range_start = calc_stopping_range(E_start)
        range_end = calc_stopping_range(E_end)
        slice_thickness_g_cm2 = range_start - range_end
        print(f"Slice Thickness = Range({E_start} MeV) - Range({E_end} MeV) = {range_start:.5f} - {range_end:.5f} = {slice_thickness_g_cm2:.5f} g/cm^2")
        
        # Calculate average cross section for the slice
        avg_sigma_mb = (sigma_mb[E_start] + sigma_mb[E_end]) / 2.0
        avg_sigma_cm2 = avg_sigma_mb * 1e-27 # convert millibarns to cm^2
        print(f"Average Cross-Section = ({sigma_mb[E_start]} + {sigma_mb[E_end]})/2 = {avg_sigma_mb:.2f} mb = {avg_sigma_cm2:.4e} cm^2")
        
        # Calculate production rate for this slice
        slice_production_rate = N_atoms_per_gram * slice_thickness_g_cm2 * avg_sigma_cm2 * I_protons_per_sec
        print(f"Slice Production Rate = {N_atoms_per_gram:.4e} * {slice_thickness_g_cm2:.5f} * {avg_sigma_cm2:.4e} * {I_protons_per_sec:.4e} = {slice_production_rate:.4e} atoms/s")
        total_production_rate += slice_production_rate
    
    print(f"\nTotal Production Rate (R) = {total_production_rate:.4e} atoms/s (decays/s at saturation)")
    
    # --- 4. Calculate Final Activity ---
    print("\nStep 4: Calculate Final Activity at End-of-Bombardment (EOB)")
    
    # Convert times to seconds
    t_irradiation_sec = t_irradiation_hr * 3600
    T_half_sec = T_half_days * 24 * 3600
    
    # Calculate decay constant (lambda)
    lambda_per_sec = math.log(2) / T_half_sec
    print(f"Decay Constant (λ) = ln(2) / ({T_half_days} * 24 * 3600 s) = {lambda_per_sec:.4e} s^-1")
    
    # Calculate saturation factor
    saturation_factor = 1 - math.exp(-lambda_per_sec * t_irradiation_sec)
    print(f"Saturation Factor = 1 - e^(-{lambda_per_sec:.4e} * {t_irradiation_sec}) = {saturation_factor:.5f}")
    
    # Calculate activity in Bq
    activity_bq = total_production_rate * saturation_factor
    print(f"Activity (Bq) = R * (Saturation Factor) = {total_production_rate:.4e} * {saturation_factor:.5f} = {activity_bq:.4e} Bq")

    # --- 5. Convert to millicuries ---
    print("\nStep 5: Convert Activity to Millicuries (mCi)")
    # 1 Ci = 3.7e10 Bq, so 1 mCi = 3.7e7 Bq
    activity_mci = activity_bq / 3.7e7
    print(f"Final Activity (mCi) = {activity_bq:.4e} Bq / 3.7e7 Bq/mCi = {activity_mci:.2f} mCi")
    
    # Return the final value in the required format
    return activity_mci

if __name__ == '__main__':
    final_yield = calculate_tb155_yield()
    print(f"\n<<< {final_yield:.2f} >>>")
