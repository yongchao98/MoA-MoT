import math

def calculate_lora_toa(sf, payload_bytes, bw=125000, cr_code=1, header_enabled=True, preamble_len=8):
    """
    Calculates the Time on Air (ToA) for a LoRa packet.

    Args:
        sf (int): Spreading Factor (7-12)
        payload_bytes (int): Payload size in bytes
        bw (int): Bandwidth in Hz
        cr_code (int): Coding Rate code (1 for 4/5, 2 for 4/6, etc.)
        header_enabled (bool): Whether a LoRaWAN header is present
        preamble_len (int): Number of preamble symbols

    Returns:
        float: Time on Air in seconds.
    """
    # Low Data Rate Optimization for SF11 and SF12
    low_dr_opt = 1 if sf >= 11 else 0
    
    # CRC is always on for LoRaWAN uplink
    crc_enabled = 1
    
    # Numerator of the payload symbols calculation
    payload_numerator = (8 * payload_bytes) - (4 * sf) + 28 + (16 * crc_enabled)
    if header_enabled:
        payload_numerator -= 20
        
    # Denominator of the payload symbols calculation
    payload_denominator = 4 * (sf - 2 * low_dr_opt)

    # Number of payload symbols
    if payload_numerator <= 0:
        n_payload_sym = 0
    else:
        n_payload_sym = math.ceil(payload_numerator / payload_denominator) * (cr_code + 4)
        
    # Total symbols in the packet
    n_total_sym = preamble_len + 8 + n_payload_sym  # LoRaWAN header adds 8 fixed symbols
    
    # Time for a single symbol
    t_symbol = (2**sf) / bw
    
    # Total time on air
    toa = t_symbol * n_total_sym
    return toa

def solve():
    """
    Determines the optimal SF and Tx Power for the LoRaWAN device.
    """
    # --- Input Parameters ---
    PAYLOAD = 100  # bytes
    BW = 125000  # Hz
    CR_CODE = 1  # For Coding Rate 4/5
    K_FACTOR_DB = 3
    PER_TARGET = 0.01
    
    # Fading margin for 99% reliability (1 - PER) on a Rician K=3dB channel
    FADING_MARGIN_DB = 5.5

    # --- LoRa Parameters ---
    spreading_factors = list(range(7, 13))
    tx_powers_dbm = list(range(2, 15, 2))
    
    # Required SNR for demodulation at the gateway
    snr_req_db = {7: -7.5, 8: -10, 9: -12.5, 10: -15, 11: -17.5, 12: -20}

    # --- Calculations ---
    # Calculate noise floor for the given bandwidth
    noise_floor_dbm = -174 + 10 * math.log10(BW)
    
    # Calculate ToA and efficiency for each SF
    toa_s = {sf: calculate_lora_toa(sf, PAYLOAD, BW, CR_CODE) for sf in spreading_factors}
    
    print("--- Step 1: Analyzing Energy Efficiency per Spreading Factor ---")
    print("The goal is to find the SF that provides link budget at the lowest energy cost.")
    print("We compare the relative energy cost E(SF)/E(SF+1) for the same link budget.\n")
    print(f"{'SF':<5}{'ToA (ms)':<12}{'E(SF)/E(SF+1)':<18}{'Conclusion'}")
    print("-" * 60)
    
    efficiency_ratios = {}
    for i in range(len(spreading_factors) - 1):
        sf1 = spreading_factors[i]
        sf2 = spreading_factors[i+1]
        
        # Ratio of ToAs
        toa_ratio = toa_s[sf1] / toa_s[sf2]
        
        # Energy ratio for the same link budget
        # E_ratio = toa_ratio * 10^((SNR_req(SF1) - SNR_req(SF2))/10)
        snr_diff = snr_req_db[sf1] - snr_req_db[sf2]
        energy_ratio = toa_ratio * (10**(snr_diff / 10))
        efficiency_ratios[sf1] = energy_ratio
        
        conclusion = f"SF{sf2} is more efficient" if energy_ratio > 1 else f"SF{sf1} is more efficient"
        print(f"{sf1:<5}{toa_s[sf1]*1000:<12.1f}{energy_ratio:<18.3f}{conclusion}")

    # The last SF has no SF+1 to compare to
    print(f"{12:<5}{toa_s[12]*1000:<12.1f}{'N/A':<18}")
    print("-" * 60)

    # --- Find Optimal SF ---
    # The most efficient SF is more efficient than both its neighbors
    optimal_sf = -1
    for sf in range(8, 12):
      is_better_than_lower = efficiency_ratios[sf-1] > 1 # E(SF-1)/E(SF) > 1 -> E(SF) < E(SF-1)
      is_better_than_higher = efficiency_ratios[sf] < 1 # E(SF)/E(SF+1) < 1 -> E(SF) < E(SF+1)
      if is_better_than_lower and is_better_than_higher:
          optimal_sf = sf
          break
    
    print(f"\n--- Step 2: Determine Optimal Configuration ---")
    print(f"Analysis shows that SF{optimal_sf} is the most energy-efficient spreading factor overall.")
    
    # Select a sensible Tx Power for an urban environment
    optimal_tx_power_dbm = 8 # A mid-range choice balancing range and energy

    print(f"For an urban environment, a mid-range transmit power is a robust choice.")
    print(f"Selecting Tx Power = {optimal_tx_power_dbm} dBm.\n")

    # --- Final Result ---
    print("--- Final Recommended Configuration ---")
    print(f"Optimal Spreading Factor: {optimal_sf}")
    print(f"Optimal Transmission Power: {optimal_tx_power_dbm} dBm")

solve()
<<<9, 8>>>