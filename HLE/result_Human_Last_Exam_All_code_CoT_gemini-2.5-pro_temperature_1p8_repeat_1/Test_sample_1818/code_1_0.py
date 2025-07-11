import math

def calculate_lorawan_toa(sf, pl, bw, cr_code_val, n_preamble=8):
    """Calculates the Time on Air (ToA) for a LoRaWAN packet."""
    # Low Data Rate Optimization for SF11 and SF12
    de = 0
    if sf >= 11:
        de = 1
    
    # Symbol Time
    t_sym = (2**sf) / bw

    # Preamble duration
    t_preamble = (n_preamble + 4.25) * t_sym

    # Number of payload symbols calculation
    # Formula components: Payload Length (PL), Spreading Factor (SF), Header (H=0), CRC (16), Coding Rate (CR)
    # The term (cr_code_val + 4) represents the coding rate factor (e.g., for CR=4/5, it's 5)
    payload_numerator = 8 * pl - 4 * sf + 28 + 16 
    payload_denominator = 4 * (sf - 2 * de)
    
    # Check for non-positive denominator, which can happen for SF<=6 with DE on
    if payload_denominator <= 0:
      return float('inf')

    payload_symb_nb = 8 + max(0, math.ceil(payload_numerator / payload_denominator) * cr_code_val)
    
    # Payload duration
    t_payload = payload_symb_nb * t_sym
    
    # Total ToA
    toa = t_preamble + t_payload
    return toa

def solve_lorawan_optimization():
    """
    Determines the optimal SF and TxP for energy efficiency
    based on the problem's constraints.
    """
    # Given Parameters
    PAYLOAD_BYTES = 100
    BANDWIDTH_HZ = 125000
    CODING_RATE_STR = "4/5"
    TX_POWERS_DBM = list(range(2, 15, 2))
    SPREADING_FACTORS = list(range(7, 13))

    # LoRaWAN specific parameters for calculation
    # Coding Rate 4/5 -> CR code = 1, effective multiplier = 5
    CR_CODE = 1 
    CR_MULTIPLIER = CR_CODE + 4
    
    # Assumptions for SNR Calculation
    # Baseline SNR sensitivity for 10% PER at 125kHz
    SNR_SENSITIVITY_10_PER = {
        7: -7.5, 8: -10.0, 9: -12.5, 10: -15.0, 11: -17.5, 12: -20.0
    }
    # Margin to improve PER from 10% to 1%
    PER_MARGIN_DB = 2.0
    # Fading margin for Rician channel (K=3dB) for 99% reliability
    FADING_MARGIN_DB = 6.0
    
    print("--- LoRaWAN Energy Optimization Analysis ---\n")
    print("This script determines the optimal Spreading Factor (SF) and Transmit Power (TxP)")
    print("to minimize energy while maintaining a Packet Error Rate (PER) <= 1%.")
    print("\nKey Assumptions:")
    print(f"- PER Target: 1% requires a {PER_MARGIN_DB} dB margin over standard sensitivity.")
    print(f"- Rician Fading (K=3dB): Requires a {FADING_MARGIN_DB} dB fading margin for 99% reliability.")
    print("\n--- Calculating Optimal Configuration ---\n")
    
    results = []
    
    for sf in SPREADING_FACTORS:
        # 1. Calculate Time on Air (ToA)
        toa_s = calculate_lorawan_toa(sf, PAYLOAD_BYTES, BANDWIDTH_HZ, CR_MULTIPLIER)
        
        # 2. Determine Required SNR
        snr_base = SNR_SENSITIVITY_10_PER[sf]
        required_snr_db = snr_base + PER_MARGIN_DB + FADING_MARGIN_DB
        
        # 3. Calculate Energy Factor (proportional to energy)
        # Energy ~ P_tx * ToA. P_tx is proportional to 10^(SNR/10).
        # So, Energy Factor ~ 10^(SNR/10) * ToA
        power_factor = 10**(required_snr_db / 10)
        energy_factor = power_factor * toa_s
        
        results.append({
            "sf": sf,
            "toa_ms": toa_s * 1000,
            "required_snr_db": required_snr_db,
            "energy_factor": energy_factor
        })

    # Print the results in a formatted table
    print(f"{'SF':<5} | {'ToA (ms)':<12} | {'Req. SNR (dB)':<15} | {'Relative Energy Factor':<25}")
    print("-" * 65)
    for res in results:
        print(f"SF{res['sf']:<2} | {res['toa_ms']:<12.1f} | {res['required_snr_db']:<15.1f} | {res['energy_factor']:<25.2f}")

    # 4. Find the optimal SF with the minimum energy factor
    optimal_config = min(results, key=lambda x: x['energy_factor'])
    optimal_sf = optimal_config['sf']
    
    # 5. Determine the optimal TxP
    # ADR will choose the lowest power that meets the link budget.
    # Without a specified path loss, we select the lowest available power
    # as representative of an optimal (good link) scenario.
    optimal_tx_power_dbm = TX_POWERS_DBM[0]

    print("\n--- Conclusion ---\n")
    print(f"The most energy-efficient Spreading Factor is SF{optimal_sf}, as it has the lowest Relative Energy Factor.")
    print(f"The Adaptive Data Rate (ADR) mechanism aims to use the minimum power necessary.")
    print(f"For a device with a sufficient link budget, the optimal transmission power would be the lowest possible setting.")
    print(f"\nOptimal Configuration: Spreading Factor = {optimal_sf}, Transmission Power = {optimal_tx_power_dbm} dBm")
    
    # 6. Calculate the final energy value for the optimal configuration
    optimal_toa_s = optimal_config['toa_ms'] / 1000.0
    optimal_tx_power_watts = (10**((optimal_tx_power_dbm - 30) / 10))
    energy_joules = optimal_tx_power_watts * optimal_toa_s

    print("\nThe energy consumption for this optimal transmission is calculated as:")
    print("Energy (J) = Transmit Power (W) * Time on Air (s)")
    print(f"Energy = {optimal_tx_power_watts:.4f} W * {optimal_toa_s:.4f} s = {energy_joules:.6f} J ({energy_joules*1000:.3f} mJ)")

solve_lorawan_optimization()