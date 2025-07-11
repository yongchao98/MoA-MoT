import math

def calculate_toa_and_energy():
    """
    Calculates the optimal SF and TP for a LoRaWAN device based on energy efficiency
    for a given link budget requirement.
    """
    # --- 1. Define Parameters ---
    PL = 100  # Payload in bytes
    BW = 125000  # Bandwidth in Hz
    CR_code = 1  # Coding rate 4/5 corresponds to CR_code=1 in formula
    N_preamble = 8  # Number of preamble symbols
    H = 1  # Header is enabled

    # Available Spreading Factors and Transmit Powers
    spreading_factors = list(range(7, 13))
    transmit_powers_dbm = list(range(2, 15, 2))

    # SNR requirements and fading margin
    snr_req_awgn = {7: -7.5, 8: -10.0, 9: -12.5, 10: -15.0, 11: -17.5, 12: -20.0}
    fading_margin = 4.5  # For Rician K=3dB at 99% reliability
    
    # Noise calculation
    noise_floor_dbm = -174 + 10 * math.log10(BW)
    gateway_noise_figure = 3.0
    total_noise = noise_floor_dbm + gateway_noise_figure
    
    # Assumed required link budget for an urban environment
    required_link_budget = 136.0

    # --- 2. Analysis ---
    results = []
    for sf in spreading_factors:
        # Calculate Time on Air (ToA)
        DE = 1 if sf >= 11 else 0  # Data rate optimization for SF11, SF12
        t_symbol = (2**sf) / BW
        
        # Number of payload symbols calculation
        payload_numerator = 8 * PL - 4 * sf + 28 + 16 - 20 * H
        payload_denominator = 4 * (sf - 2 * DE)
        n_payload = 8 + max(math.ceil(payload_numerator / payload_denominator) * (CR_code + 4), 0)
        
        toa_s = (N_preamble + 4.25 + n_payload) * t_symbol

        # Calculate required SNR including fading margin
        snr_req_fading = snr_req_awgn[sf] + fading_margin
        
        # Calculate gateway sensitivity
        gateway_sensitivity = snr_req_fading + total_noise

        # Find minimum TP to meet link budget
        min_tp_needed = required_link_budget + gateway_sensitivity
        
        min_tp_available = None
        for tp in transmit_powers_dbm:
            if tp >= min_tp_needed:
                min_tp_available = tp
                break
        
        # If a valid TP is found, calculate energy and store results
        if min_tp_available is not None:
            power_mw = 10**((min_tp_available - 30) / 10) * 1000
            energy_mj = power_mw * toa_s
            actual_link_budget = min_tp_available - gateway_sensitivity
            results.append({
                "sf": sf,
                "tp_dbm": min_tp_available,
                "toa_ms": toa_s * 1000,
                "energy_mj": energy_mj,
                "link_budget_db": actual_link_budget
            })

    # --- 3. Output Results ---
    print(f"Analysis for a required Link Budget of {required_link_budget} dB")
    print("-" * 70)
    print(f"{'SF':<5} | {'Min TP (dBm)':<15} | {'ToA (ms)':<12} | {'Energy (mJ)':<15} | {'Achieved LB (dB)':<18}")
    print("-" * 70)

    min_energy = float('inf')
    optimal_config = None

    for res in results:
        print(f"SF{res['sf']:<4} | {res['tp_dbm']:<15.1f} | {res['toa_ms']:<12.2f} | {res['energy_mj']:<15.2f} | {res['link_budget_db']:<18.1f}")
        if res['energy_mj'] < min_energy:
            min_energy = res['energy_mj']
            optimal_config = res
            
    print("-" * 70)
    
    if optimal_config:
        print("\nOptimal Configuration Found:")
        sf = optimal_config['sf']
        tp = optimal_config['tp_dbm']
        toa_ms = optimal_config['toa_ms']
        power_mw = 10**((tp - 30) / 10) * 1000
        energy_mj = optimal_config['energy_mj']
        
        print(f"The most energy-efficient configuration to meet the requirement is SF{sf} at {tp} dBm.")
        print("\nFinal Energy Calculation:")
        print(f"  Transmit Power ({tp} dBm)      = {power_mw:.2f} mW")
        print(f"  Time on Air (SF{sf})           = {toa_ms:.2f} ms")
        print(f"  Energy = {power_mw:.2f} mW * {toa_ms/1000:.4f} s = {energy_mj:.2f} mJ")
        
        # Final answer in the required format
        print(f"\n<<<SF{sf}, {tp} dBm>>>")
    else:
        print("No configuration could meet the required link budget.")
        print("\n<<<Infeasible>>>")


calculate_toa_and_energy()