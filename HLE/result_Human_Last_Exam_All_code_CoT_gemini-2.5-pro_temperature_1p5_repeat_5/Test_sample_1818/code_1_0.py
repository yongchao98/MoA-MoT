import math

# Plan:
# 1. Define all known parameters of the LoRaWAN system and the environment.
# 2. Model the requirements. This involves determining the required Signal-to-Noise Ratio (SNR) at the receiver
#    to meet the 1% Packet Error Rate (PER) in a Rician fading channel. We will do this by taking the known
#    SNR limit for each Spreading Factor (SF) and adding a Fading Margin.
# 3. Assume a realistic Path Loss (PL) for the urban link. A crucial assumption is that the link is challenging enough
#    to require a trade-off between SF and Transmit Power (TxP). We will select a PL that makes the lowest SF (SF7)
#    viable only at the highest available power. This forces a meaningful optimization problem.
# 4. Iterate through all available Spreading Factors (SF7 to SF12).
# 5. For each SF, calculate:
#    a. The Time on Air (ToA) for the 100-byte payload.
#    b. The required average receiver SNR.
#    c. The minimum required TxP to overcome the assumed path loss and noise.
#    d. The actual TxP to use, which is the lowest available power level that meets or exceeds the requirement.
#    e. The total energy consumed for one transmission (Energy = Power * ToA).
# 6. Compare the energy consumption for all viable (SF, TxP) pairs.
# 7. Identify the pair with the minimum energy consumption as the optimal setting.

def solve_lora_optimization():
    """
    Calculates the optimal SF and TxP for a LoRaWAN device to minimize energy.
    """
    # --- Step 1: Define Parameters ---

    # LoRaWAN Parameters
    PAYLOAD_BYTES = 100
    BANDWIDTH_HZ = 125000
    CODING_RATE_DENOM = 5
    PREAMBLE_SYMBOLS = 8
    AVAILABLE_TX_POWER_DBM = list(range(2, 16, 2))  # 2 to 14 in 2dB steps
    SPREADING_FACTORS = list(range(7, 13))

    # Environmental & Link Parameters
    K_FACTOR_DB = 3.0  # Rician K-factor
    TARGET_PER = 0.01
    # Assumption: A fading margin is needed to ensure 99% reliability (1 - PER).
    # For a Rician channel with K=3dB, a margin of ~5.5dB is appropriate for 99% reliability.
    FADING_MARGIN_DB = 5.5

    # Assumption: A typical gateway receiver noise figure.
    GATEWAY_NOISE_FIGURE_DB = 6.0

    # Assumption: Path Loss is set to a value that creates a meaningful trade-off.
    # We choose a path loss where SF7 is only just viable with the maximum available power.
    PATH_LOSS_DB = 133.0

    print(f"--- Assumptions ---")
    print(f"Path Loss (PL): {PATH_LOSS_DB} dB")
    print(f"Gateway Noise Figure (NF): {GATEWAY_NOISE_FIGURE_DB} dB")
    print(f"Rician Fading Margin for 99% reliability (1% PER): {FADING_MARGIN_DB} dB\n")

    # --- Step 2: Model Requirements (SNR & Noise Floor) ---

    # SNR limits (sensitivity) for LoRa modulation in AWGN
    SNR_LIMITS_DB = {
        7: -7.5, 8: -10.0, 9: -12.5, 10: -15.0, 11: -17.5, 12: -20.0
    }

    # Required average SNR at receiver = SNR limit + Fading Margin
    required_snr_db = {sf: limit + FADING_MARGIN_DB for sf, limit in SNR_LIMITS_DB.items()}

    # Calculate receiver noise floor
    noise_floor_dbm = -174 + 10 * math.log10(BANDWIDTH_HZ) + GATEWAY_NOISE_FIGURE_DB

    print(f"--- System Calculations ---")
    print(f"Receiver Noise Floor: {noise_floor_dbm:.2f} dBm\n")


    # --- Step 3 & 4: Iterate and Calculate ---

    results = []

    def calculate_time_on_air(sf, payload_bytes):
        """Calculates the Time on Air for a given LoRa configuration."""
        low_dr_opt = 1 if sf >= 11 else 0
        cr_code = CODING_RATE_DENOM - 4
        t_symbol = (2**sf) / BANDWIDTH_HZ
        numerator = 8 * payload_bytes - 4 * sf + 28 + 16
        denominator = 4 * (sf - 2 * low_dr_opt)
        n_payload_sym = 8 + math.ceil(numerator / denominator) * (cr_code + 4)
        n_total_sym = PREAMBLE_SYMBOLS + 4.25 + n_payload_sym
        toa_s = n_total_sym * t_symbol
        return toa_s

    print("--- Evaluating each Spreading Factor ---")
    for sf in SPREADING_FACTORS:
        toa = calculate_time_on_air(sf, PAYLOAD_BYTES)
        req_snr = required_snr_db[sf]
        req_tx_power_dbm = PATH_LOSS_DB + req_snr + noise_floor_dbm
        
        chosen_tx_power_dbm = -1
        for p in AVAILABLE_TX_POWER_DBM:
            if p >= req_tx_power_dbm:
                chosen_tx_power_dbm = p
                break
        
        if chosen_tx_power_dbm == -1:
            print(f"SF{sf}: Required TxP ({req_tx_power_dbm:.2f} dBm) exceeds max available power. Not viable.")
            continue

        chosen_tx_power_watts = 10**((chosen_tx_power_dbm - 30) / 10)
        energy_joules = chosen_tx_power_watts * toa

        results.append({
            "SF": sf,
            "Chosen TxP (dBm)": chosen_tx_power_dbm,
            "ToA (s)": toa,
            "Energy (mJ)": energy_joules * 1000,
        })

        print(f"SF{sf}: Req. TxP={req_tx_power_dbm:.2f} dBm -> Use {chosen_tx_power_dbm} dBm. "
              f"ToA={toa:.3f} s. Energy={energy_joules*1000:.2f} mJ")

    # --- Step 6 & 7: Compare results and find the optimum ---
    if not results:
        print("\nNo viable SF/TxP combination found for the given path loss.")
    else:
        optimal_config = min(results, key=lambda x: x["Energy (mJ)"])

        print("\n--- Summary & Conclusion ---")
        print(f"The optimal configuration to minimize energy is:")
        final_sf = optimal_config['SF']
        final_txp = optimal_config['Chosen TxP (dBm)']
        final_energy = optimal_config['Energy (mJ)']
        final_toa = optimal_config['ToA (s)']
        
        print(f"Spreading Factor: {final_sf}")
        print(f"Transmission Power: {final_txp} dBm")
        print(f"\nThis configuration's energy consumption ({final_energy:.2f} mJ) is calculated as:")
        print(f"Energy (mJ) = (10^(({final_txp} - 30) / 10)) * {final_toa:.3f} * 1000")
        
        # Final answer in requested format
        print(f"\n<<<Optimal SF: {final_sf}, Optimal TxP: {final_txp} dBm>>>")

if __name__ == '__main__':
    solve_lora_optimization()