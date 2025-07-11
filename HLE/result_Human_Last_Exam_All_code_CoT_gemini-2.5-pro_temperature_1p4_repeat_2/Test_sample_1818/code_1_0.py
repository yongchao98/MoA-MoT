import math

def solve_lorawan_optimization():
    """
    Calculates the optimal LoRaWAN Spreading Factor and Transmission Power
    to minimize energy for a given scenario.
    """
    
    # --- 1. Define Constants and Parameters ---
    
    # LoRaWAN Parameters from the problem description
    PAYLOAD_BYTES = 100
    BANDWIDTH_HZ = 125000
    CODING_RATE_NUM = 1  # For CR = 4/5
    
    # LoRaWAN protocol constants
    PREAMBLE_SYMBOLS = 8
    HEADER_ENABLED = True  # Explicit header means H=0 in the formula
    CRC_ENABLED = True     # CRC is on
    
    # RF and Channel Parameter Assumptions (based on problem description and typical values)
    # Fading margin for 99% reliability (PER < 1%) on a Rician channel with K=3dB.
    # This is a key assumption derived from channel models.
    FADING_MARGIN_DB = 4.5
    
    # A representative path loss for an urban environment. The problem is unsolvable
    # without this assumption, as it dictates the required power.
    PATH_LOSS_DB = 135
    
    # A typical gateway receiver noise floor (includes thermal noise and receiver noise figure).
    NOISE_FLOOR_DBM = -120
    
    # Available device transmit power levels as per the problem
    TX_POWER_LEVELS_DBM = list(range(2, 16, 2))  # [2, 4, 6, 8, 10, 12, 14]
    
    # Spreading Factors and their baseline SNR requirements for demodulation in an AWGN channel
    SF_SNR_MAP = {
        7: -7.5,
        8: -10.0,
        9: -12.5,
        10: -15.0,
        11: -17.5,
        12: -20.0
    }
    
    # --- 2. Helper Functions ---
    
    def calculate_toa_seconds(sf):
        """Calculates the Time on Air (ToA) in seconds for a given SF."""
        H = 0 if HEADER_ENABLED else 1
        CRC = 1 if CRC_ENABLED else 0
        DE = 1 if sf >= 11 else 0  # Data rate optimization for SF11, SF12
        
        t_symbol = (2**sf) / BANDWIDTH_HZ
        
        # LoRaWAN specified formula for number of payload symbols
        numerator = 8 * PAYLOAD_BYTES - 4 * sf + 28 + 16 * CRC - 20 * H
        denominator = 4 * (sf - 2 * DE)
        
        payload_symbol_count = 8 + math.ceil(max(0, numerator / denominator)) * (CODING_RATE_NUM + 4)
        
        total_symbols = PREAMBLE_SYMBOLS + 4.25 + payload_symbol_count
        return total_symbols * t_symbol

    def find_actual_tx_power(required_dbm):
        """Finds the lowest available Tx power that meets the required power."""
        if required_dbm > max(TX_POWER_LEVELS_DBM):
            return None  # Requirement cannot be met
        
        for power in TX_POWER_LEVELS_DBM:
            if power >= required_dbm:
                return power
        # This fallback should not be reached due to the initial check
        return max(TX_POWER_LEVELS_DBM)

    # --- 3. Main Calculation Logic ---
    
    results = []
    
    print("--- LoRaWAN Energy Optimization Analysis ---")
    print(f"Key Assumptions:\n  - Path Loss: {PATH_LOSS_DB} dB\n  - Fading Margin for PER < 1%: {FADING_MARGIN_DB} dB\n  - Gateway Noise Floor: {NOISE_FLOOR_DBM} dBm\n")
    
    for sf, snr_awgn in SF_SNR_MAP.items():
        # a. Calculate final SNR requirement including fading margin
        snr_required_db = snr_awgn + FADING_MARGIN_DB
        
        # b. Calculate required Tx Power using the link budget equation:
        # Tx_Power = SNR_req + Path_Loss + Noise_Floor
        required_tx_power_dbm = snr_required_db + PATH_LOSS_DB + NOISE_FLOOR_DBM
        
        # c. Find the actual Tx Power to be assigned by the server
        actual_tx_power_dbm = find_actual_tx_power(required_tx_power_dbm)
        
        print(f"Analyzing SF{sf}...")
        
        if actual_tx_power_dbm is not None:
            # d. Calculate Time on Air (ToA)
            toa_s = calculate_toa_seconds(sf)
            
            # e. Calculate Energy Consumption: Energy = Power (W) * Time (s)
            power_watts = (10**(actual_tx_power_dbm / 10)) / 1000
            energy_joules = power_watts * toa_s
            
            results.append({
                "sf": sf,
                "actual_tx_power_dbm": actual_tx_power_dbm,
                "energy_joules": energy_joules
            })
            print(f"  > Valid Configuration: Power = {actual_tx_power_dbm} dBm, ToA = {toa_s * 1000:.1f} ms, Energy = {energy_joules * 1000:.1f} mJ")
        else:
            print(f"  > Invalid: Required power ({required_tx_power_dbm:.1f} dBm) exceeds maximum available power ({max(TX_POWER_LEVELS_DBM)} dBm).")

    # --- 4. Find the Optimal Setting and Print Final Answer ---
    
    if not results:
        print("\nError: No valid SF/Tx Power combination found for the assumed path loss.")
    else:
        optimal_setting = min(results, key=lambda x: x['energy_joules'])
        
        final_sf = optimal_setting['sf']
        final_tx = optimal_setting['actual_tx_power_dbm']
        
        print("\n--- Optimal Configuration ---")
        print("To minimize energy while ensuring PER < 1%, the network server should assign:")
        print(f"Final Equation:")
        print(f"Optimal Spreading Factor = {final_sf}")
        print(f"Optimal Transmission Power = {final_tx} dBm")
        print(f"Resulting Minimum Energy = {optimal_setting['energy_joules'] * 1000:.1f} mJ")

solve_lorawan_optimization()