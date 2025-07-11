import math

def calculate_lora_toa_ms(sf, payload_length, bw_khz=125, coding_rate_den=5, preamble_len=8, header_enabled=True, crc_on=True):
    """
    Calculates the Time on Air (ToA) for a LoRa packet.

    Args:
        sf (int): Spreading Factor (7-12).
        payload_length (int): Payload length in bytes.
        bw_khz (int): Bandwidth in kHz (e.g., 125, 250, 500).
        coding_rate_den (int): Denominator of the coding rate (e.g., 5 for 4/5, 8 for 4/8).
        preamble_len (int): Number of preamble symbols.
        header_enabled (bool): Whether the LoRaWAN header is present.
        crc_on (bool): Whether CRC is enabled in the payload.

    Returns:
        float: Time on Air in milliseconds.
    """
    # Symbol duration in ms
    t_sym = (2**sf) / (bw_khz)

    # Low Data Rate Optimization (LDRO) is enabled for SF > 10 on 125kHz
    de = 1 if sf > 10 and bw_khz == 125 else 0
    
    # Coding Rate value for the formula (1 for 4/5, 2 for 4/6, etc.)
    cr_val = max(coding_rate_den - 4, 1)

    # Number of preamble symbols
    t_preamble = (preamble_len + 4.25) * t_sym
    
    # Header is 20 bytes, but the formula uses a bit flag 'H'
    # H = 0 for explicit header, H = 1 for no header.
    h = 0 if header_enabled else 1
    
    # CRC presence flag
    crc = 1 if crc_on else 0

    # Numerator of the payload symbols calculation
    payload_symb_nb_num = 8 * payload_length - 4 * sf + 28 + 16 * crc - 20 * h
    
    # Denominator of the payload symbols calculation
    payload_symb_nb_den = 4 * (sf - 2 * de)
    
    if payload_symb_nb_den == 0:
        return float('inf') # Avoid division by zero
    
    # Number of payload symbols
    payload_symb_nb = 8 + max(math.ceil(payload_symb_nb_num / payload_symb_nb_den) * cr_val, 0)

    t_payload = payload_symb_nb * t_sym

    return t_preamble + t_payload

def solve_lorawan_energy():
    """
    Solves for the optimal SF and Tx Power to minimize energy consumption.
    """
    # --- Input Parameters ---
    PAYLOAD_BYTES = 100
    BANDWIDTH_KHZ = 125
    CODING_RATE_DEN = 5 # for a 4/5 coding rate
    
    # --- System Parameters ---
    AVAILABLE_TX_POWERS_DBM = [2, 4, 6, 8, 10, 12, 14]
    SPREADING_FACTORS = [7, 8, 9, 10, 11, 12]
    
    # SNR sensitivity (dB) for PER < 1% at 125kHz BW
    SNR_SENSITIVITY_DB = [-7.5, -10, -12.5, -15, -17.5, -20]
    
    # --- Channel and Link Budget Assumptions ---
    # Rician K-factor of 3dB requires a fading margin for 99% reliability.
    # An engineering approximation of 4.5 dB is used.
    RICIAN_FADING_MARGIN_DB = 4.5
    
    # Assume a typical urban path loss and gateway noise floor.
    # The absolute values are less important than their difference between SFs.
    PATH_LOSS_DB = 130
    NOISE_FLOOR_DBM = -120 # For 125kHz BW including gateway noise figure

    # --- Calculation ---
    results = []
    min_energy_mj = float('inf')
    optimal_sf = None
    optimal_tx_power_dbm = None
    optimal_toa_s = None
    optimal_power_mw = None

    print("--- Analyzing Energy Consumption for Each Spreading Factor ---")
    
    for i, sf in enumerate(SPREADING_FACTORS):
        # 1. Calculate Time on Air (ToA)
        toa_ms = calculate_lora_toa_ms(sf, PAYLOAD_BYTES, bw_khz=BANDWIDTH_KHZ, coding_rate_den=CODING_RATE_DEN)
        toa_s = toa_ms / 1000.0

        # 2. Determine Required Transmit Power
        base_snr_req_db = SNR_SENSITIVITY_DB[i]
        total_snr_req_db = base_snr_req_db + RICIAN_FADING_MARGIN_DB
        
        required_tx_power_dbm = total_snr_req_db - NOISE_FLOOR_DBM + PATH_LOSS_DB
        
        # 3. Select lowest available Tx Power that meets the requirement
        actual_tx_power_dbm = None
        for pwr in AVAILABLE_TX_POWERS_DBM:
            if pwr >= required_tx_power_dbm:
                actual_tx_power_dbm = pwr
                break
        
        # If no available power level is sufficient, this SF is not viable
        if actual_tx_power_dbm is None:
            print(f"SF{sf:<2}: Not viable. Required TxPower ({required_tx_power_dbm:.1f} dBm) exceeds max available power.")
            continue
            
        # 4. Calculate Energy Consumption
        tx_power_mw = 10**(actual_tx_power_dbm / 10)
        energy_mj = tx_power_mw * toa_s
        
        results.append({
            "sf": sf,
            "toa_s": toa_s,
            "req_tx_dbm": required_tx_power_dbm,
            "actual_tx_dbm": actual_tx_power_dbm,
            "energy_mj": energy_mj
        })

        print(f"SF{sf:<2}: ToA={toa_s:>6.3f}s, Req TxPwr={required_tx_power_dbm:>5.1f}dBm -> Use {actual_tx_power_dbm}dBm. Energy: {energy_mj:.3f} mJ")
        
        # 5. Check for new minimum
        if energy_mj < min_energy_mj:
            min_energy_mj = energy_mj
            optimal_sf = sf
            optimal_tx_power_dbm = actual_tx_power_dbm
            optimal_toa_s = toa_s
            optimal_power_mw = tx_power_mw

    # --- Final Result ---
    print("\n--- Optimal Configuration ---")
    if optimal_sf is not None:
        print(f"The optimal combination to minimize energy is Spreading Factor {optimal_sf} at {optimal_tx_power_dbm} dBm.")
        print(f"This configuration achieves the required Packet Error Rate with the lowest energy consumption.")
        print("\nFinal Energy Calculation:")
        # The final required format: output each number in the final equation
        print(f"{min_energy_mj:.3f} mJ = {optimal_power_mw:.2f} mW * {optimal_toa_s:.3f} s")
    else:
        print("No viable configuration found for the given parameters.")

    return optimal_sf, optimal_tx_power_dbm

# Execute the function to find the solution
solve_lorawan_energy()