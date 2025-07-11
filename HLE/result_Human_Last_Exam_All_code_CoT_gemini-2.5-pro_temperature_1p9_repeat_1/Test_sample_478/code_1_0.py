import math

def solve_acoustics_problem():
    """
    This script calculates the Signal-to-Noise Ratio (SNR) based on the provided acoustic data.
    It follows these main steps:
    1.  Uses mixed sound level measurements to determine the source strength of the train and construction.
    2.  Calculates the noise level from both sources at the location of the people.
    3.  Combines these noise levels to find a total noise figure.
    4.  Calculates the SNR by subtracting the total noise level from the signal level.
    """

    # --- Step 1: Define Given Values ---
    # Signal level at the source (people talking at their location, interpreted as r=1m)
    L_signal = 75.0  # dB

    # Measurement Point 1 (Train + People)
    L_combo1 = 100.0         # Total dB
    dist_train_to_p1 = 10.0  # meters
    dist_people_to_p1 = 20.0 # meters

    # Measurement Point 2 (Construction + People)
    L_combo2 = 115.0         # Total dB
    dist_constr_to_p2 = 20.0 # meters
    dist_people_to_p2 = 30.0 # meters

    # Location of Interest (where people are talking)
    dist_train_to_people = 30.0   # meters
    dist_constr_to_people = 50.0  # meters
    
    # --- Step 2: Calculate Source Strengths (at a reference distance of 1m) ---

    # 2a. Find the source strength of the train
    L_people_at_p1 = L_signal - 20 * math.log10(dist_people_to_p1)
    I_combo1 = 10**(L_combo1 / 10.0)
    I_people_at_p1 = 10**(L_people_at_p1 / 10.0)
    I_train_at_p1 = I_combo1 - I_people_at_p1
    L_train_at_p1 = 10 * math.log10(I_train_at_p1)
    L_train_1m = L_train_at_p1 + 20 * math.log10(dist_train_to_p1)

    # 2b. Find the source strength of the construction
    L_people_at_p2 = L_signal - 20 * math.log10(dist_people_to_p2)
    I_combo2 = 10**(L_combo2 / 10.0)
    I_people_at_p2 = 10**(L_people_at_p2 / 10.0)
    I_constr_at_p2 = I_combo2 - I_people_at_p2
    L_constr_at_p2 = 10 * math.log10(I_constr_at_p2)
    L_constr_1m = L_constr_at_p2 + 20 * math.log10(dist_constr_to_p2)

    # --- Step 3: Calculate Total Noise Level at the People's Location ---

    # 3a. Calculate individual noise levels at the people's location
    L_train_at_people = L_train_1m - 20 * math.log10(dist_train_to_people)
    L_constr_at_people = L_constr_1m - 20 * math.log10(dist_constr_to_people)
    
    # 3b. Combine noise levels
    I_train_noise = 10**(L_train_at_people / 10.0)
    I_constr_noise = 10**(L_constr_at_people / 10.0)
    I_total_noise = I_train_noise + I_constr_noise
    L_total_noise = 10 * math.log10(I_total_noise)

    # --- Step 4: Calculate SNR and Print the Final Equation ---
    snr_db = L_signal - L_total_noise

    print("The final calculation is SNR = S - N, where S is the signal level and N is the total noise level.")
    print("\n--- Final Equation Components ---")
    print(f"Signal Level (S): {L_signal:.2f} dB")
    print(f"Total Noise Level (N): {L_total_noise:.2f} dB")
    print(f"Final Equation: SNR = {L_signal:.2f} - {L_total_noise:.2f}")
    print(f"\nResulting Signal-to-Noise Ratio (SNR): {snr_db:.2f} dB")
    
# Run the calculation
solve_acoustics_problem()