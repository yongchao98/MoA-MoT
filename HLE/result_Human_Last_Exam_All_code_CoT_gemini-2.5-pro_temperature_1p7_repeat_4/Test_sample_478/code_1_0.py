import math

def solve_snr_problem():
    """
    Calculates the Signal-to-Noise Ratio (SNR) based on the provided acoustic data.
    """
    
    # --- Helper functions for acoustic calculations ---
    
    def db_to_intensity(L):
        """Converts sound level in dB to a normalized linear intensity."""
        return 10**(L / 10.0)

    def intensity_to_db(I):
        """Converts a normalized linear intensity back to sound level in dB."""
        if I <= 0:
            return -float('inf') # Intensity must be positive for log
        return 10.0 * math.log10(I)

    def get_level_at_distance(L_ref, r_ref, r_new):
        """Calculates sound level at a new distance from a source using inverse square law."""
        return L_ref - 20.0 * math.log10(r_new / r_ref)

    # --- Problem Parameters ---
    
    # Signal level at the measurement location
    L_signal = 75.0  # dB

    # Assumed source strength of people talking (at 1 meter)
    L_people_source_1m = 75.0  # dB

    # Measurement Point 1: Train + People
    L_total_1 = 100.0  # dB
    r_train_1 = 10.0   # meters
    r_people_1 = 20.0  # meters

    # Measurement Point 2: Construction + People
    L_total_2 = 115.0  # dB
    r_construction_2 = 20.0 # meters
    r_people_2 = 30.0   # meters
    
    # Final location (where people/listener are)
    r_train_final = 30.0      # meters
    r_construction_final = 50.0 # meters
    
    # --- Step 1: Find Train's Source Strength ---
    
    # Calculate the sound level of the people at measurement point 1
    L_people_at_1 = get_level_at_distance(L_people_source_1m, 1.0, r_people_1)
    
    # Convert total level and people's level to intensity to subtract them
    I_total_1 = db_to_intensity(L_total_1)
    I_people_at_1 = db_to_intensity(L_people_at_1)
    
    # The people's sound intensity is negligible compared to the total, but we calculate precisely.
    I_train_at_1 = I_total_1 - I_people_at_1
    L_train_at_10m = intensity_to_db(I_train_at_1)
    
    # --- Step 2: Find Construction's Source Strength ---

    # Calculate the sound level of the people at measurement point 2
    L_people_at_2 = get_level_at_distance(L_people_source_1m, 1.0, r_people_2)
    
    # Convert total level and people's level to intensity
    I_total_2 = db_to_intensity(L_total_2)
    I_people_at_2 = db_to_intensity(L_people_at_2)
    
    # The people's sound is again negligible, but we calculate precisely.
    I_construction_at_2 = I_total_2 - I_people_at_2
    L_construction_at_20m = intensity_to_db(I_construction_at_2)
    
    # --- Step 3: Calculate Total Noise at the Final Location ---
    
    # Calculate the noise level from the train at the listener's location
    L_noise_train_final = get_level_at_distance(L_train_at_10m, r_train_1, r_train_final)
    
    # Calculate the noise level from construction at the listener's location
    L_noise_construction_final = get_level_at_distance(L_construction_at_20m, r_construction_2, r_construction_final)

    # Convert individual noise levels to intensities, add them, and convert back to a total noise level
    I_noise_train_final = db_to_intensity(L_noise_train_final)
    I_noise_construction_final = db_to_intensity(L_noise_construction_final)
    
    I_noise_total_final = I_noise_train_final + I_noise_construction_final
    L_noise_total_final = intensity_to_db(I_noise_total_final)
    
    # --- Step 4: Calculate the final SNR ---
    
    SNR = L_signal - L_noise_total_final
    
    # --- Print Results ---
    
    print("Step 1: The calculated noise level from the train at the final location is {:.2f} dB.".format(L_noise_train_final))
    print("Step 2: The calculated noise level from the construction site at the final location is {:.2f} dB.".format(L_noise_construction_final))
    print("Step 3: The total combined noise level is {:.2f} dB.".format(L_noise_total_final))
    print("\nFinal Calculation:")
    print("The Signal-to-Noise Ratio (SNR) is the Signal Level minus the Total Noise Level.")
    # The final equation with each number printed out
    print("SNR = {:.2f} dB - {:.2f} dB = {:.2f} dB".format(L_signal, L_noise_total_final, SNR))

solve_snr_problem()