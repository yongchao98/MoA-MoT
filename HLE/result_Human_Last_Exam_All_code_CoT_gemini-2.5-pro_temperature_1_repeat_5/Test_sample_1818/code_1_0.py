import math

def solve_lorawan_optimization():
    """
    Calculates the optimal LoRaWAN Spreading Factor and Transmit Power
    to minimize energy consumption while meeting a PER requirement.
    """
    # --- 1. Define Constants and Parameters ---
    PAYLOAD_BYTES = 100
    BANDWIDTH = 125000  # Hz
    CODING_RATE_INFO = {'name': '4/5', 'code': 1}
    PREAMBLE_SYMBOLS = 8
    IMPLICIT_HEADER = False # Using an explicit header
    CRC_ENABLED = True

    AVAILABLE_TX_POWERS_DBM = sorted([2, 4, 6, 8, 10, 12, 14])
    SPREADING_FACTORS = range(7, 13)  # SF7 to SF12

    # --- 2. Environment and Link Model Parameters ---
    # For urban Rician fading with K=3dB, a fading margin is needed.
    # A 5 dB margin is a reasonable choice for a Packet Error Rate of ~1%.
    RICIAN_FADING_MARGIN_DB = 5.0

    # Base SNR requirements (demodulation threshold in an AWGN channel)
    SNR_REQUIREMENTS_AWGN = {
        7: -7.5, 8: -10.0, 9: -12.5, 10: -15.0, 11: -17.5, 12: -20.0
    }

    # A representative path loss for an urban environment.
    PATH_LOSS_DB = 135.0

    # --- Helper Functions ---
    def calculate_noise_floor_dbm(bw):
        # Thermal noise power in dBm: -174 dBm/Hz + 10*log10(Bandwidth)
        return -174 + 10 * math.log10(bw)

    def calculate_toa_ms(sf, payload, bw, cr_code, preamble, implicit_header, crc):
        # Symbol duration in ms
        t_symbol = (2**sf / bw) * 1000
        
        # Preamble duration
        t_preamble = (preamble + 4.25) * t_symbol

        # Payload calculation terms
        ih = 1 if implicit_header else 0
        de = 1 if sf in [11, 12] and bw == 125000 else 0
        crc_val = 1 if crc else 0
        cr_multiplier = cr_code + 4

        # Number of payload symbols
        numerator = 8 * payload - 4 * sf + 28 + 16 * crc_val - 20 * ih
        denominator = 4 * (sf - 2 * de)
        
        num_payload_symbols = 8 + max(math.ceil(numerator / denominator) * cr_multiplier, 0)
        
        # Payload duration
        t_payload = num_payload_symbols * t_symbol
        
        return t_preamble + t_payload

    # --- Main Calculation Logic ---
    print("Analyzing LoRaWAN Parameters for Optimal Energy Consumption\n")
    print(f"Assumptions:\n- Path Loss: {PATH_LOSS_DB} dB\n- Rician Fading Margin for 1% PER: {RICIAN_FADING_MARGIN_DB} dB\n")

    results = []
    noise_floor = calculate_noise_floor_dbm(BANDWIDTH)

    for sf in SPREADING_FACTORS:
        # Calculate the total required SNR at the receiver
        required_snr = SNR_REQUIREMENTS_AWGN[sf] + RICIAN_FADING_MARGIN_DB
        
        # Calculate the required transmit power using the link budget equation
        # TxPower = RequiredSNR - RxGain + PathLoss + NoiseFloor + ImplementationMargin - TxGain
        # Simplified: TxPower_req = RequiredSNR + PathLoss - NoiseFloor
        required_tx_power_dbm = required_snr + PATH_LOSS_DB + noise_floor

        # Find the lowest available Tx power that meets the requirement
        chosen_tx_power_dbm = None
        for power in AVAILABLE_TX_POWERS_DBM:
            if power >= required_tx_power_dbm:
                chosen_tx_power_dbm = power
                break
        
        if chosen_tx_power_dbm is None:
            # This SF is not viable with available power levels
            continue

        # Calculate Time on Air for this SF
        toa_ms = calculate_toa_ms(
            sf=sf, payload=PAYLOAD_BYTES, bw=BANDWIDTH,
            cr_code=CODING_RATE_INFO['code'], preamble=PREAMBLE_SYMBOLS,
            implicit_header=IMPLICIT_HEADER, crc=CRC_ENABLED
        )

        # Calculate relative energy consumption (proportional to mJ)
        # Energy is proportional to Power (mW) * Time (s)
        power_mw = 10**(chosen_tx_power_dbm / 10)
        relative_energy = power_mw * toa_ms
        
        results.append({
            'sf': sf,
            'required_snr': required_snr,
            'required_tx_power': required_tx_power_dbm,
            'chosen_tx_power': chosen_tx_power_dbm,
            'toa_ms': toa_ms,
            'relative_energy': relative_energy
        })

    # --- Display Results and Find Optimum ---
    print("--- Step-by-Step Analysis ---")
    for res in results:
        print(f"\nSF{res['sf']}:")
        print(f"  - Required SNR: {SNR_REQUIREMENTS_AWGN[res['sf']]:.1f} (base) + {RICIAN_FADING_MARGIN_DB:.1f} (margin) = {res['required_snr']:.1f} dB")
        print(f"  - Required Tx Power: {res['required_snr']:.1f} + {PATH_LOSS_DB:.1f} - {abs(noise_floor):.1f} = {res['required_tx_power']:.1f} dBm")
        print(f"  - Chosen Tx Power: {res['chosen_tx_power']} dBm (next available level)")
        print(f"  - Time on Air: {res['toa_ms']:.2f} ms")
        print(f"  - Relative Energy Equation: 10^({res['chosen_tx_power']}/10) * {res['toa_ms']:.2f} = {res['relative_energy']:.1f}")

    if not results:
        print("\nCould not find a viable configuration with the given parameters.")
        return

    optimal_choice = min(results, key=lambda x: x['relative_energy'])
    
    print("\n--- Conclusion ---")
    print("The optimal configuration is the one with the lowest relative energy score.")
    print(f"The best choice is SF{optimal_choice['sf']} with a transmit power of {optimal_choice['chosen_tx_power']} dBm.")

    global final_answer
    final_answer = f"SF{optimal_choice['sf']}, {optimal_choice['chosen_tx_power']} dBm"

# Execute the analysis
solve_lorawan_optimization()
# The final answer is captured in a global variable to be formatted at the end.
print(f"\n<<<SF{optimal_choice['sf']}, {optimal_choice['chosen_tx_power']} dBm>>>")