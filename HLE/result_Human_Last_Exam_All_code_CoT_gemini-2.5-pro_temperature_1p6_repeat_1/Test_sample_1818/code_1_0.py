import math

def calculate_energy_consumption():
    """
    Calculates the optimal LoRaWAN SF and Tx power for minimal energy consumption
    in a specified urban environment.
    """

    # 1. Define LoRaWAN Parameters
    PL = 100  # Payload in bytes
    BW = 125000  # Bandwidth in Hz
    CR_CODE = 1  # Coding Rate 4/5 corresponds to code 1
    N_PREAMBLE = 8  # Number of preamble symbols
    CRC = True  # CRC is on
    IH = False # Implicit header is off for SF7-12
    
    # Available operating parameters
    spreading_factors = list(range(7, 13))
    tx_powers_dbm = list(range(2, 15, 2))

    # Receiver sensitivity (required SNR for PER < 1%) for each SF
    snr_req = {
        7: -7.5,
        8: -10,
        9: -12.5,
        10: -15,
        11: -17.5,
        12: -20,
    }

    # 2. Model the Communication Link
    PATH_LOSS_DB = 140  # Assumed path loss for a challenging urban environment
    NOISE_FIGURE_DB = 6   # Typical gateway receiver noise figure
    # Noise Floor (dBm) = -174 + 10*log10(BW) + Noise Figure
    noise_floor_dbm = -174 + 10 * math.log10(BW) + NOISE_FIGURE_DB

    print(f"System Parameters:")
    print(f"  - Payload: {PL} bytes")
    print(f"  - Bandwidth: {BW/1000} kHz")
    print(f"  - Coding Rate: 4/5")
    print(f"  - Assumed Urban Path Loss: {PATH_LOSS_DB} dB")
    print(f"  - Receiver Noise Floor: {noise_floor_dbm:.2f} dBm\n")

    # Store results to find the minimum
    results = []
    
    # 3. & 4. Iterate through each SF to find the optimal configuration
    print("Analyzing Energy Consumption for each Spreading Factor...")
    for sf in spreading_factors:
        # Calculate Time on Air (ToA)
        # Low data rate optimization is enabled for SF11 and SF12 on 125kHz
        de = 1 if sf >= 11 else 0
        t_sym = (2**sf) / BW
        
        # Preamble duration
        t_preamble = (N_PREAMBLE + 4.25) * t_sym
        
        # Payload duration
        payload_part_numerator = 8 * PL - 4 * sf + 28 + (16 if CRC else 0) - (20 if IH else 0)
        payload_part_denominator = 4 * (sf - 2 * de)
        payload_symb_nb = 8 + max(0, math.ceil(payload_part_numerator / payload_part_denominator) * (CR_CODE + 4))
        
        t_payload = payload_symb_nb * t_sym
        toa_s = t_preamble + t_payload

        # Calculate required transmit power
        # Tx_Power(dBm) >= Path_Loss(dB) + Noise_Floor(dBm) + SNR_req(dB)
        req_tx_power_dbm = PATH_LOSS_DB + noise_floor_dbm + snr_req[sf]

        # Check if this SF is viable with available power levels
        if req_tx_power_dbm > max(tx_powers_dbm):
            print(f"\nSF{sf}:")
            print(f"  - Required Tx Power ({req_tx_power_dbm:.2f} dBm) exceeds maximum available power ({max(tx_powers_dbm)} dBm).")
            print(f"  - Verdict: NOT VIABLE for this path loss.")
            continue

        # Find the next available discrete power level
        actual_tx_power_dbm = 0
        for p in tx_powers_dbm:
            if p >= req_tx_power_dbm:
                actual_tx_power_dbm = p
                break
        
        # Calculate energy consumption
        actual_tx_power_mw = 10**(actual_tx_power_dbm / 10)
        energy_mj = actual_tx_power_mw * toa_s * 1000 # To get mJ

        results.append({
            "sf": sf,
            "toa_ms": toa_s * 1000,
            "req_tx_dbm": req_tx_power_dbm,
            "actual_tx_dbm": actual_tx_power_dbm,
            "actual_tx_mw": actual_tx_power_mw,
            "energy_mj": energy_mj
        })
        
        print(f"\nSF{sf}:")
        print(f"  - Time on Air: {toa_s*1000:.2f} ms")
        print(f"  - Required Tx Power: {req_tx_power_dbm:.2f} dBm")
        print(f"  - Assigned Tx Power: {actual_tx_power_dbm} dBm ({actual_tx_power_mw:.2f} mW)")
        print(f"  - Energy Consumed = {toa_s*1000:.2f} ms * {actual_tx_power_mw:.2f} mW = {energy_mj:.2f} mJ")

    # 5. Determine the Winner
    if not results:
        print("\nNo viable configuration found for the given path loss.")
        return
        
    optimal_config = min(results, key=lambda x: x['energy_mj'])

    print("\n------------------------------------------------------------")
    print("                     CONCLUSION")
    print("------------------------------------------------------------")
    print("To achieve the PER requirement with the lowest energy consumption for the given urban environment,")
    print("the network server should assign the following parameters:")
    
    sf_opt = int(optimal_config['sf'])
    tx_opt_dbm = int(optimal_config['actual_tx_dbm'])
    
    print("\nOptimal Configuration:")
    print(f"  Spreading Factor: SF{sf_opt}")
    print(f"  Transmission Power: {tx_opt_dbm} dBm")

    # Final equation format as requested
    print("\nFinal Energy Calculation for Optimal Setting:")
    print(f"  Energy = Time on Air * Transmit Power")
    print(f"  Energy = {optimal_config['toa_ms']:.2f} ms * {optimal_config['actual_tx_mw']:.2f} mW = {optimal_config['energy_mj']:.2f} mJ")
    
    # Required final answer format
    print(f"\n<<<Spreading Factor: SF{sf_opt}, Transmission Power: {tx_opt_dbm} dBm>>>")


if __name__ == "__main__":
    calculate_energy_consumption()
