import math
from collections import Counter

def solve_lorawan_optimization():
    """
    Determines the optimal Spreading Factor and Transmit Power for a LoRaWAN device
    to minimize energy while maintaining a PER < 1%.
    """
    
    # --- 1. System and Channel Parameters ---
    PAYLOAD_BYTES = 100
    BANDWIDTH_HZ = 125000
    CODING_RATE_CODE = 1  # Represents CR 4/5
    
    # Environment and Reliability
    PER_TARGET = 0.01
    # For a Rician channel with K=3dB, a fading margin of ~4.5 dB is required for 99% reliability.
    FADING_MARGIN_DB = 4.5
    
    # Receiver and Channel Noise
    THERMAL_NOISE_PER_HZ_DBM = -174
    NOISE_FLOOR_DBM = THERMAL_NOISE_PER_HZ_DBM + 10 * math.log10(BANDWIDTH_HZ)

    # Device Options
    SPREADING_FACTORS = list(range(7, 13))
    TX_POWERS_DBM = list(range(2, 15, 2))
    
    # Required SNR for demodulation at the receiver for each SF
    SNR_REQ_DB = {7: -7.5, 8: -10, 9: -12.5, 10: -15, 11: -17.5, 12: -20}

    # --- 2. Time on Air (ToA) Calculation ---
    def calculate_toa(sf):
        t_symbol = (2**sf) / BANDWIDTH_HZ
        ldro_enabled = t_symbol > 0.016  # LDRO is on for symbol times > 16ms (SF11, SF12 at 125kHz)
        de = 1 if ldro_enabled else 0
        
        # Payload symbols calculation based on Semtech datasheets
        numerator = 8 * PAYLOAD_BYTES - 4 * sf + 28 + 16
        denominator = 4 * (sf - 2 * de)
        n_payload_sym = 8 + max(0, math.ceil(numerator / denominator) * (CODING_RATE_CODE + 4))
        
        # Total ToA = Preamble Time + Payload Time
        n_preamble = 8
        total_symbols = n_preamble + 4.25 + n_payload_sym
        toa_seconds = total_symbols * t_symbol
        return toa_seconds

    toas = {sf: calculate_toa(sf) for sf in SPREADING_FACTORS}

    # --- 3. Link Budget Requirement Calculation ---
    min_rssi_req = {
        sf: SNR_REQ_DB[sf] + FADING_MARGIN_DB + NOISE_FLOOR_DBM
        for sf in SPREADING_FACTORS
    }
    
    # --- 4. ADR Simulation across Path Losses ---
    optimal_choices = []
    # Simulate for path losses from 125 dB to 153 dB in 0.1 dB steps
    path_loss_range_db = [i * 0.1 for i in range(1250, 1531)] 

    for path_loss_db in path_loss_range_db:
        min_energy_mj = float('inf')
        best_setting_for_pl = None

        for sf in SPREADING_FACTORS:
            # Calculate the transmit power required to overcome the path loss for this SF
            required_tx_power_dbm = path_loss_db + min_rssi_req[sf]

            # Find the lowest available transmit power that meets the requirement
            chosen_tx_power_dbm = -1
            for pwr in TX_POWERS_DBM:
                if pwr >= required_tx_power_dbm:
                    chosen_tx_power_dbm = pwr
                    break
            
            # If a valid power level exists for this SF
            if chosen_tx_power_dbm != -1:
                tx_power_mw = 10**(chosen_tx_power_dbm / 10)
                energy_mj = tx_power_mw * toas[sf]

                if energy_mj < min_energy_mj:
                    min_energy_mj = energy_mj
                    best_setting_for_pl = (sf, chosen_tx_power_dbm)
        
        if best_setting_for_pl:
            optimal_choices.append(best_setting_for_pl)

    # --- 5. Identify the Most Optimal Setting ---
    if not optimal_choices:
        print("No optimal setting found within the specified path loss range.")
        return

    # Count the occurrences of each winning setting
    winning_counts = Counter(optimal_choices)
    # Find the setting that won the most times
    overall_winner, count = winning_counts.most_common(1)[0]
    
    winner_sf, winner_pwr_dbm = overall_winner
    winner_toa = toas[winner_sf]
    winner_pwr_mw = 10**(winner_pwr_dbm / 10)
    winner_energy_mj = winner_pwr_mw * winner_toa

    print("--- LoRaWAN Energy Optimization Result ---")
    print(f"\nThe optimal setting that minimizes energy over the widest range of conditions is:")
    print(f"\nSpreading Factor: SF{winner_sf}")
    print(f"Transmission Power: {winner_pwr_dbm} dBm\n")

    print("--- Final Energy Calculation for this Setting ---")
    print(f"Time on Air for SF{winner_sf} (100-byte payload): {winner_toa:.4f} seconds")
    print(f"Transmission Power ({winner_pwr_dbm} dBm): {winner_pwr_mw:.2f} mW")
    print("\nFinal Equation:")
    print(f"Energy = Power * Time on Air")
    print(f"Energy = {winner_pwr_mw:.2f} mW * {winner_toa:.4f} s = {winner_energy_mj:.2f} mJ")

solve_lorawan_optimization()