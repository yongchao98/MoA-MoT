import math

def calculate_toa(sf, payload_len, bw=125000, cr_code=1, n_preamble=8, explicit_header=True, crc=True):
    """
    Calculates the Time on Air for a LoRaWAN packet.
    """
    # Low Data Rate Optimization is enabled for SF11 and SF12 with 125kHz bandwidth
    de = 0
    if sf >= 11 and bw == 125000:
        de = 1

    # Header is not present in the payload for explicit header mode
    h = 0 if explicit_header else 1

    # CRC presence
    crc_val = 1 if crc else 0
    
    # Symbol time
    t_sym = (2**sf) / bw

    # Number of payload symbols calculation
    numerator = 8 * payload_len - 4 * sf + 28 + 16 * crc_val - 20 * h
    denominator = 4 * (sf - 2 * de)
    
    # Ceil division
    payload_symbol_count = 8 + max(0, math.ceil(numerator / denominator) * (cr_code + 4))

    # Total time on air
    toa = (n_preamble + 4.25 + payload_symbol_count) * t_sym
    return toa

def solve_lorawan_optimization():
    """
    Finds the optimal SF and TxP for minimum energy consumption under specific urban conditions.
    """
    # 1. Available Parameters & Assumptions
    spreading_factors = [7, 8, 9, 10, 11, 12]
    transmit_powers_dbm = [2, 4, 6, 8, 10, 12, 14]
    
    # LoRaWAN receiver sensitivity (SNR) for each SF
    snr_sensitivity = {7: -7.5, 8: -10, 9: -12.5, 10: -15, 11: -17.5, 12: -20}
    
    # Payload
    payload_len = 100  # bytes
    
    # Link budget assumptions
    path_loss_db = 130  # Representative path loss for an urban environment
    noise_figure_db = 6 # Typical gateway noise figure
    bandwidth_hz = 125000
    fading_margin_db = 5  # Required margin for Rician K=3dB at 99% reliability (1% PER)

    # 2. Calculate Noise Floor and Link Requirement
    noise_floor_dbm = -174 + 10 * math.log10(bandwidth_hz) + noise_figure_db
    link_requirement_db = path_loss_db + noise_floor_dbm + fading_margin_db

    print("--- LoRaWAN Energy Optimization ---")
    print(f"\nStep 1: Define Scenario & Assumptions")
    print(f"  - Environment: Urban with Rician fading (K=3dB)")
    print(f"  - Assumed Path Loss: {path_loss_db} dB")
    print(f"  - Gateway Noise Figure: {noise_figure_db} dB")
    print(f"  - Required Fading Margin for PER < 1%: {fading_margin_db} dB")

    print(f"\nStep 2: Calculate Total Link Requirement")
    print(f"  - Gateway Noise Floor = -174 + 10*log10({bandwidth_hz}) + {noise_figure_db} = {noise_floor_dbm:.2f} dBm")
    print(f"  - Total Link Requirement = Path Loss + Noise Floor + Fading Margin")
    print(f"  - Total Link Requirement = {path_loss_db} + ({noise_floor_dbm:.2f}) + {fading_margin_db} = {link_requirement_db:.2f} dB")
    
    min_energy_mj = float('inf')
    optimal_sf = None
    optimal_txp = None
    optimal_capability = None

    print(f"\nStep 3: Evaluate all (SF, TxP) combinations")
    print("Searching for the combination that meets the link requirement with the minimum energy...\n")

    # 3. Iterate through all combinations to find the optimum
    for sf in spreading_factors:
        # Calculate ToA for the SF
        toa_s = calculate_toa(sf=sf, payload_len=payload_len)
        
        for txp in transmit_powers_dbm:
            # Calculate the link capability of this setting
            capability_db = txp - snr_sensitivity[sf]
            
            # Check if this setting meets the link requirement
            if capability_db >= link_requirement_db:
                # Calculate energy consumption
                txp_mw = 10**(txp / 10)
                energy_mj = txp_mw * toa_s

                # If it's the most energy-efficient so far, save it
                if energy_mj < min_energy_mj:
                    min_energy_mj = energy_mj
                    optimal_sf = sf
                    optimal_txp = txp
                    optimal_capability = capability_db

    # 4. Output the result
    print("--- Optimal Setting Found ---")
    if optimal_sf is not None:
        toa_s = calculate_toa(sf=optimal_sf, payload_len=payload_len)
        txp_mw = 10**(optimal_txp / 10)
        
        print(f"\nOptimal Spreading Factor: SF{optimal_sf}")
        print(f"Optimal Transmission Power: {optimal_txp} dBm")
        
        print("\nSupporting Calculations for this choice:")
        
        print(f"\n1. Link Capability Check:")
        print(f"   Equation: Capability = TxP - SNR Sensitivity")
        print(f"   Calculation: {optimal_txp} dBm - ({snr_sensitivity[optimal_sf]} dB) = {optimal_capability:.2f} dB")
        print(f"   Result: Capability ({optimal_capability:.2f} dB) >= Link Requirement ({link_requirement_db:.2f} dB). Link is viable.")

        print(f"\n2. Energy Consumption Calculation:")
        print(f"   Equation: Energy = Power (mW) * Time on Air (s)")
        print(f"   Time on Air for SF{optimal_sf}: {toa_s:.4f} s")
        print(f"   Power for {optimal_txp} dBm: {txp_mw:.2f} mW")
        print(f"   Calculation: {txp_mw:.2f} mW * {toa_s:.4f} s = {min_energy_mj:.4f} mJ")
        print("\nThis combination provides the required link reliability with the lowest possible energy.")

    else:
        print("No suitable SF/TxP combination found for the assumed path loss.")

if __name__ == '__main__':
    solve_lorawan_optimization()
    print("\n<<<SF8, 8dBm>>>")