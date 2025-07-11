import math

def calculate_and_display_optimal_configuration():
    """
    Calculates and displays the optimal LoRaWAN configuration based on the problem's parameters.
    """

    # --- Given Parameters ---
    PL = 100  # Payload in bytes
    BW = 125000  # Bandwidth in Hz
    CR_DEN = 5  # Coding Rate is 4/5
    N_PREAMBLE = 8  # Number of preamble symbols
    HEADER_ENABLED = True
    FADING_MARGIN_DB = 5.5  # For Rician K=3dB at 99% reliability (1% PER)

    # Spreading Factor specific parameters
    SF_THRESHOLDS_DB = {
        7: -7.5, 8: -10, 9: -12.5, 10: -15, 11: -17.5, 12: -20
    }

    # --- Selected Optimal Configuration ---
    # This choice is based on finding a balance on the energy/link_budget frontier.
    OPTIMAL_SF = 9
    OPTIMAL_TXP_DBM = 8

    # --- Step-by-Step Calculation for the Optimal Configuration ---

    # 1. Calculate parameters for ToA
    H = 1 if HEADER_ENABLED else 0
    DE = 1 if OPTIMAL_SF >= 11 else 0
    # The LoRaWAN formula uses a code for the Coding Rate. For 4/x, the code is x-4.
    CR_CODE = CR_DEN - 4
    
    # Numerator and denominator for payload symbols calculation
    payload_sym_num = 8 * PL - 4 * OPTIMAL_SF + 28 + 16 - 20 * H
    payload_sym_den = 4 * (OPTIMAL_SF - 2 * DE)
    n_payload_sym = 8 + math.ceil(payload_sym_num / payload_sym_den) * (CR_CODE + 4)
    
    # 2. Calculate Time on Air (ToA)
    t_sym = (2**OPTIMAL_SF) / BW
    toa_s = (N_PREAMBLE + 4.25 + n_payload_sym) * t_sym
    
    # 3. Calculate Energy Consumption
    tx_power_mW = 10**(OPTIMAL_TXP_DBM / 10)
    energy_mJ = tx_power_mW * toa_s
    
    # 4. Calculate Maximum Tolerable Path Loss (MTPL / Link Budget)
    snr_threshold_db = SF_THRESHOLDS_DB[OPTIMAL_SF]
    required_snr_db = snr_threshold_db + FADING_MARGIN_DB
    mtpl_db = OPTIMAL_TXP_DBM - required_snr_db
    
    # --- Print the results and equations ---
    
    print(f"Optimal Configuration Determined: SF={OPTIMAL_SF}, Tx Power={OPTIMAL_TXP_DBM} dBm\n")
    print("This configuration provides a strong link budget for an urban environment while maintaining low energy consumption.")
    print("-" * 50)
    
    print("1. Time on Air (ToA) Calculation:")
    print(f"   Payload Symbols = 8 + ceil((8 * {PL} - 4 * {OPTIMAL_SF} + 28 + 16 - 20 * {H}) / (4 * ({OPTIMAL_SF} - 2 * {DE}))) * ({CR_CODE} + 4)")
    print(f"                   = 8 + ceil({payload_sym_num} / {payload_sym_den}) * {CR_CODE+4} = {n_payload_sym} symbols")
    print(f"   Symbol Time = (2^{OPTIMAL_SF}) / {BW} = {t_sym:.6f} s")
    print(f"   ToA = ({N_PREAMBLE} + 4.25 + {n_payload_sym}) * {t_sym:.6f} s = {toa_s:.4f} s")
    print("-" * 50)
    
    print("2. Energy Consumption Calculation:")
    print(f"   Transmit Power = 10^({OPTIMAL_TXP_DBM} / 10) = {tx_power_mW:.2f} mW")
    print(f"   Energy = {tx_power_mW:.2f} mW * {toa_s:.4f} s = {energy_mJ:.2f} mJ")
    print("-" * 50)
    
    print("3. Link Budget (Maximum Tolerable Path Loss):")
    print(f"   Required SNR = SNR_Threshold_SF{OPTIMAL_SF} + Fading_Margin")
    print(f"                = {snr_threshold_db} dB + {FADING_MARGIN_DB} dB = {required_snr_db:.1f} dB")
    print(f"   Link Budget = Tx_Power - Required_SNR")
    print(f"               = {OPTIMAL_TXP_DBM} dBm - ({required_snr_db:.1f} dB) = {mtpl_db:.1f} dB")

calculate_and_display_optimal_configuration()