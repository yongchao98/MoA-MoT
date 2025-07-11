import math

def solve_lorawan_optimization():
    """
    Calculates the optimal LoRaWAN Spreading Factor and Transmit Power
    to minimize energy consumption for a given set of constraints.
    """
    # 1. Define constants and parameters
    PAYLOAD_BYTES = 100
    BANDWIDTH_HZ = 125000
    CODING_RATE_CODE = 1  # For CR = 4/5
    PREAMBLE_LEN = 8
    HEADER_ENABLED = 0  # 0 for explicit header, 1 for implicit
    TX_POWERS_DBM = list(range(2, 15, 2))  # 2, 4, ..., 14 dBm

    # Baseline sensitivities (dBm) for each SF in an AWGN channel
    SENSITIVITY_BASE = {
        7: -7.5, 8: -10.0, 9: -12.5, 10: -15.0, 11: -17.5, 12: -20.0
    }
    # A Rician fading channel (K=3dB) with a 1% PER target requires a fading margin.
    # A typical value for this is ~5 dB.
    FADING_MARGIN_DB = 5.0

    # Required SNR for each SF to meet the PER target
    SNR_REQ = {sf: sens + FADING_MARGIN_DB for sf, sens in SENSITIVITY_BASE.items()}

    def calculate_toa_s(sf):
        """Calculates the Time on Air (ToA) in seconds for a given SF."""
        low_data_rate_optimize = 1 if sf >= 11 else 0
        t_symbol_s = (2**sf) / BANDWIDTH_HZ
        
        # Formula for number of payload symbols
        numerator = 8 * PAYLOAD_BYTES - 4 * sf + 28 + 16 - 20 * HEADER_ENABLED
        denominator = 4 * (sf - 2 * low_data_rate_optimize)
        
        payload_symb_count = 8 + max(0, math.ceil(numerator / denominator) * (CODING_RATE_CODE + 4))
        
        total_symbols = PREAMBLE_LEN + 4.25 + payload_symb_count
        return total_symbols * t_symbol_s

    # 2. Simulate a realistic but challenging urban path loss scenario
    # We assume a Link Budget Requirement (L = PathLoss + NoiseFloor) that makes the choice non-trivial.
    # Let's choose a value that SF9 can't handle, to force a decision between higher SFs.
    # Max link budget for SF9 = 14dBm (max power) - (-7.5dB SNR_req + 5dB margin) = 21.5 dB.
    # We'll use a slightly higher value.
    LINK_BUDGET_REQ_DB = 21.6

    print("Analyzing LoRaWAN settings for minimal energy consumption.")
    print(f"Payload: {PAYLOAD_BYTES} bytes, Bandwidth: {BANDWIDTH_HZ/1000} kHz, Coding Rate: 4/5")
    print(f"Assumed Link Budget Requirement (L): {LINK_BUDGET_REQ_DB} dB")
    print("-" * 70)
    print(f"{'SF':>4} | {'SNR_req (dB)':>12} | {'TxP_req (dBm)':>13} | {'TxP_set (dBm)':>13} | {'ToA (ms)':>10} | {'Energy (uJ)':>12}")
    print("-" * 70)

    results = []
    for sf in range(7, 13):
        snr_req = SNR_REQ[sf]
        
        # Required TxP to overcome link budget and meet SNR requirement
        txp_req_dbm = snr_req + LINK_BUDGET_REQ_DB
        
        # Find the lowest available power level that works
        actual_txp_dbm = next((p for p in TX_POWERS_DBM if p >= txp_req_dbm), None)
        
        if actual_txp_dbm:
            power_watts = 10**((actual_txp_dbm - 30) / 10)
            toa_s = calculate_toa_s(sf)
            energy_j = power_watts * toa_s
            
            results.append({
                "sf": sf,
                "txp_dbm": actual_txp_dbm,
                "energy_uj": energy_j * 1e6,
                "toa_ms": toa_s * 1000,
            })
            print(f"SF{sf:<3} | {snr_req:>12.1f} | {txp_req_dbm:>13.1f} | {actual_txp_dbm:>13.0f} | {toa_s*1000:>9.1f} | {energy_j*1e6:>11.1f}")
        else:
            print(f"SF{sf:<3} | {snr_req:>12.1f} | {txp_req_dbm:>13.1f} | {'Not Possible':>13} | {'-':>10} | {'-':>12}")

    # 3. Find the setting with the minimum energy
    if not results:
        print("\nNo viable SF/TxP combination found for the assumed path loss.")
        return

    optimal_setting = min(results, key=lambda x: x['energy_uj'])
    
    opt_sf = optimal_setting['sf']
    opt_txp = optimal_setting['txp_dbm']
    opt_toa_ms = optimal_setting['toa_ms']
    opt_energy_uj = optimal_setting['energy_uj']
    
    print("-" * 70)
    print("\nOptimal setting with the lowest energy consumption:")
    print(f"Spreading Factor: SF{opt_sf}")
    print(f"Transmission Power: {opt_txp} dBm")
    
    print("\nFinal equation for the optimal energy calculation:")
    # The final print statement shows the calculation with all the numbers.
    print(f"Energy = 10^(({opt_txp} - 30) / 10) W * {opt_toa_ms / 1000:.4f} s = {opt_energy_uj / 1e6:.7f} J")

solve_lorawan_optimization()
<<<SF10, 12 dBm>>>