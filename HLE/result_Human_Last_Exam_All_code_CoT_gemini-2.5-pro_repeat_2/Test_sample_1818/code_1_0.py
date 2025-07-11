import math

def calculate_toa(sf, payload_bytes, bw=125000, cr_val=1, n_preamble=8, explicit_header=True):
    """
    Calculates the Time on Air (ToA) for a LoRa packet.
    sf: Spreading Factor (7-12)
    payload_bytes: Payload size in bytes
    bw: Bandwidth in Hz
    cr_val: Coding rate value (1 for 4/5, 2 for 4/6, etc.)
    """
    # Low Data Rate Optimization (LDRO) is used for SF11 and SF12 on 125kHz bandwidth
    ldro = 1 if sf >= 11 else 0
    
    # Symbol time in seconds
    t_symbol = (2**sf) / bw
    
    # Preamble duration
    t_preamble = (n_preamble + 4.25) * t_symbol
    
    # Header configuration
    header_on = 1 if explicit_header else 0

    # Payload symbol calculation
    numerator = 8 * payload_bytes - 4 * sf + 28 + 16 - 20 * (1 - header_on)
    denominator = 4 * (sf - 2 * ldro)
    
    # Ensure numerator is positive before ceiling operation
    payload_symbol_count_part = 0
    if numerator > 0:
        payload_symbol_count_part = math.ceil(numerator / denominator) * (cr_val + 4)
    
    n_payload_symbols = 8 + max(0, payload_symbol_count_part)
    
    # Payload duration
    t_payload = n_payload_symbols * t_symbol
    
    # Total ToA in milliseconds
    total_toa_s = t_preamble + t_payload
    return total_toa_s * 1000

def dbm_to_mw(dbm):
    """Converts dBm to milliwatts."""
    return 10**(dbm / 10)

def solve_lorawan_optimization():
    """
    Solves the LoRaWAN optimization problem to find the most energy-efficient
    SF and Tx Power combination for a given scenario.
    """
    # --- Given Parameters ---
    PAYLOAD_BYTES = 100
    FADING_MARGIN_DB = 6.4  # For Rician K=3dB at 99% reliability (1% PER)
    CODING_RATE_VAL = 1 # for 4/5
    
    # Available Spreading Factors and their minimum SNR sensitivity (dB)
    SF_SENSITIVITY = {
        7: -7.5,
        8: -10.0,
        9: -12.5,
        10: -15.0,
        11: -17.5,
        12: -20.0,
    }

    # Available Transmit Power levels (dBm)
    TX_POWERS_DBM = [2, 4, 6, 8, 10, 12, 14]

    # --- Step 1: Assume a representative Path Loss Budget to overcome ---
    # This represents the total link loss the configuration must overcome.
    # We choose a value that allows multiple SF/Tx combinations to be valid.
    PATH_LOSS_BUDGET_DB = 19.0
    
    print(f"Objective: Find the lowest energy configuration for a path loss budget of {PATH_LOSS_BUDGET_DB} dB")
    print("-" * 70)
    print(f"{'SF':<5} {'Req. SNR':<10} {'ToA (ms)':<12} {'Req. Tx (dBm)':<15} {'Actual Tx (dBm)':<17} {'Energy (mJ)':<15}")
    print("-" * 70)

    optimal_config = None
    min_energy = float('inf')

    # --- Main Calculation Loop ---
    for sf in sorted(SF_SENSITIVITY.keys()):
        # Calculate required SNR for 1% PER
        snr_required = SF_SENSITIVITY[sf] + FADING_MARGIN_DB
        
        # Calculate Time on Air for 100-byte payload
        toa_ms = calculate_toa(sf, PAYLOAD_BYTES, cr_val=CODING_RATE_VAL)

        # Determine the minimum Tx Power needed to meet the path loss budget
        # Link Budget = Tx_Power - SNR_Required >= Path_Loss_Budget
        # => Tx_Power >= Path_Loss_Budget + SNR_Required
        min_required_tx_dbm = PATH_LOSS_BUDGET_DB + snr_required
        
        # Find the smallest available Tx Power that meets the requirement
        actual_tx_dbm = None
        for p in TX_POWERS_DBM:
            if p >= min_required_tx_dbm:
                actual_tx_dbm = p
                break
        
        # If a valid Tx power is found, calculate the energy
        if actual_tx_dbm is not None:
            tx_power_mw = dbm_to_mw(actual_tx_dbm)
            energy_mj = tx_power_mw * (toa_ms / 1000)
            
            print(f"SF{sf:<2} {snr_required:<10.1f} {toa_ms:<12.3f} {min_required_tx_dbm:<15.1f} {actual_tx_dbm:<17} {energy_mj:<15.3f}")

            # Check if this configuration is the new optimum
            if energy_mj < min_energy:
                min_energy = energy_mj
                optimal_config = {
                    "sf": sf,
                    "tx_power_dbm": actual_tx_dbm,
                    "tx_power_mw": tx_power_mw,
                    "toa_ms": toa_ms,
                    "energy_mj": energy_mj
                }
        else:
            print(f"SF{sf:<2} {snr_required:<10.1f} {toa_ms:<12.3f} {min_required_tx_dbm:<15.1f} {'Not Possible':<17} {'N/A':<15}")

    # --- Print the final result ---
    print("-" * 70)
    if optimal_config:
        print("\nOptimal Configuration Found:")
        print(f"Spreading Factor: SF{optimal_config['sf']}")
        print(f"Transmission Power: {optimal_config['tx_power_dbm']} dBm")
        print("\nThis configuration provides the required link reliability with the minimum energy consumption.")
        print("\nFinal Energy Calculation:")
        print(f"Energy = Power * Time on Air")
        final_eq = (f"Energy = {optimal_config['tx_power_mw']:.2f} mW * {optimal_config['toa_ms']:.3f} ms = "
                    f"{optimal_config['energy_mj']:.3f} mJ")
        print(final_eq)
        # Final answer format for the system
        print(f"\n<<<SF{optimal_config['sf']} at {optimal_config['tx_power_dbm']} dBm>>>")
    else:
        print("No possible configuration found for the given path loss budget.")

solve_lorawan_optimization()