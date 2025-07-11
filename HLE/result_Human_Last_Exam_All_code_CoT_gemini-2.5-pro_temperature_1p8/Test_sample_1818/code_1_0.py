import math
import scipy.stats

# --- 1. LoRaWAN and Channel Parameters ---
PAYLOAD_BYTES = 100
BANDWIDTH_HZ = 125000
CODING_RATE_CODE = 1 # 1 for CR=4/5
N_PREAMBLE = 8
HEADER_ON = 0 # 0 for explicit header
CRC_ON = 1
K_FACTOR_DB = 3
TARGET_PER = 0.01

# --- 2. Device and Network Parameters ---
TX_POWERS_DBM = list(range(2, 15, 2))
SPREADING_FACTORS = list(range(7, 13))

# Baseline SNR required for demodulation per SF (from datasheets)
SNR_REQ_DB = {
    7: -7.5, 8: -10.0, 9: -12.5, 10: -15.0, 11: -17.5, 12: -20.0
}

# --- 3. Assumptions for Calculation ---
# A typical noise figure for a commercial LoRa gateway
GATEWAY_NOISE_FIGURE_DB = 6.0
# Assume a representative path loss for an urban environment
ASSUMED_PATH_LOSS_DB = 127.0


def calculate_toa_ms(sf):
    """Calculates the Time on Air (ToA) in milliseconds for a given SF."""
    # Low data rate optimization for SF11 and SF12
    de = 0
    if sf >= 11:
        de = 1
        
    t_symbol_ms = (2**sf / BANDWIDTH_HZ) * 1000
    t_preamble_ms = (N_PREAMBLE + 4.25) * t_symbol_ms
    
    # This formula is derived from the Semtech LoRa documentation
    numerator = (8 * PAYLOAD_BYTES - 4 * sf + 28 + 16 * CRC_ON - 20 * HEADER_ON)
    denominator = 4 * (sf - 2 * de)
    
    # Note: CR_CODE = 1 for 4/5, 2 for 4/6 etc. Multiplier is (CR_CODE + 4)
    payload_symbol_count = 8 + max(0, math.ceil(numerator / denominator) * (CODING_RATE_CODE + 4))
    
    t_payload_ms = payload_symbol_count * t_symbol_ms
    toa_ms = t_preamble_ms + t_payload_ms
    return toa_ms

def calculate_fading_margin_db(per_target, k_factor_db):
    """
    Calculates the required fading margin to achieve a PER target on a Rician channel.
    This uses the inverse cumulative distribution function (CDF) of a non-central
    chi-squared distribution, which models Rician fading SNR.
    """
    k_linear = 10**(k_factor_db / 10)
    # The degrees of freedom for LoRa signal envelope
    df = 2
    # The non-centrality parameter of the distribution
    nc = 2 * k_linear
    
    # We need to find the required average SNR (snr_avg) such that the
    # probability of the instantaneous SNR (snr_inst) dropping below the
    # demodulation threshold (snr_req) is equal to our target PER.
    # PER = P(snr_inst < snr_req)
    # The ratio of average to required SNR gives the margin.
    
    # ppf is the Percent Point Function (inverse of CDF)
    # The result gives a ratio that relates average SNR to the SNR at the desired percentile
    ratio_snr_req_to_avg = scipy.stats.ncx2.ppf(per_target, df, nc) / (2 * (k_linear + 1))
    
    # Margin in dB is -10*log10 of that ratio
    margin_db = -10 * math.log10(ratio_snr_req_to_avg)
    return margin_db

def find_optimal_settings():
    """Finds the optimal SF and TP for the given path loss and PER."""

    print("--- Starting Analysis ---")
    
    # Calculate fixed environmental values
    noise_floor_dbm = -174 + 10 * math.log10(BANDWIDTH_HZ) + GATEWAY_NOISE_FIGURE_DB
    fading_margin_db = calculate_fading_margin_db(TARGET_PER, K_FACTOR_DB)
    
    print(f"Assumed Urban Path Loss: {ASSUMED_PATH_LOSS_DB:.1f} dB")
    print(f"Gateway Noise Floor: {noise_floor_dbm:.1f} dBm")
    print(f"Calculated Fading Margin for PER < {TARGET_PER*100}%: {fading_margin_db:.2f} dB")
    print("-" * 25)

    min_energy_mws = float('inf')
    optimal_setting = None
    
    # Iterate through all SFs to find the best option
    for sf in SPREADING_FACTORS:
        toa_ms = calculate_toa_ms(sf)
        
        # Calculate the total required average SNR at the receiver
        snr_demod_req_db = SNR_REQ_DB[sf]
        avg_snr_req_db = snr_demod_req_db + fading_margin_db

        # Calculate the minimum TP needed to close the link with this SF
        # MinTP >= PathLoss + NoiseFloor + RequiredAvgSNR
        min_tp_req_dbm = ASSUMED_PATH_LOSS_DB + noise_floor_dbm + avg_snr_req_db
        
        # Find the first available TP that meets the requirement
        valid_tp_dbm = None
        for tp in TX_POWERS_DBM:
            if tp >= min_tp_req_dbm:
                valid_tp_dbm = tp
                break
        
        # If a valid power level exists for this SF
        if valid_tp_dbm is not None:
            power_mw = 10**(valid_tp_dbm / 10)
            energy_mws = power_mw * toa_ms
            
            print(f"Evaluating SF{sf}: Requires min TP of {min_tp_req_dbm:.2f} dBm. "
                  f"Using {valid_tp_dbm} dBm. --> Energy: {energy_mws:.0f} mWs")

            if energy_mws < min_energy_mws:
                min_energy_mws = energy_mws
                optimal_setting = {
                    "sf": sf,
                    "tp_dbm": valid_tp_dbm,
                    "toa_ms": toa_ms,
                    "power_mw": power_mw,
                    "energy_mws": energy_mws
                }
        else:
             print(f"Evaluating SF{sf}: Requires min TP of {min_tp_req_dbm:.2f} dBm. "
                   f"No available power level can meet this. Link fails.")


    # Print the final result
    print("\n--- Optimal Configuration Found ---")
    if optimal_setting:
        print(f"Optimal Spreading Factor: SF{optimal_setting['sf']}")
        print(f"Optimal Transmit Power: {optimal_setting['tp_dbm']} dBm")
        print(f"Lowest Energy Consumption: {optimal_setting['energy_mws']:.0f} mWs (milliwatt-seconds)")
        print("\nFinal Energy Calculation:")
        print(f"Energy ({optimal_setting['energy_mws']:.0f} mWs) = "
              f"Power ({optimal_setting['power_mw']:.2f} mW) * "
              f"Time on Air ({optimal_setting['toa_ms']:.2f} ms)")
    else:
        print("No viable settings found for the assumed path loss.")

if __name__ == '__main__':
    find_optimal_settings()
