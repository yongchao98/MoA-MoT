import math

def calculate_toa(payload, sf, bw=125000, cr_val=1, n_preamble=8, header=0):
    """Calculates the Time on Air for a LoRa packet."""
    # Symbol time
    t_sym = (2**sf) / bw
    
    # Low Data Rate Optimization
    de = 0
    if t_sym > 0.016:
        de = 1
        
    # Payload symbols calculation
    numerator = 8 * payload - 4 * sf + 28 + 16 - 20 * header
    denominator = 4 * (sf - 2 * de)
    pl_sym_payload = math.ceil(max(0, numerator) / denominator) * (cr_val + 4)
    
    # Total payload symbols
    pl_sym = 8 + pl_sym_payload
    
    # Preamble + Header time
    t_preamble = (n_preamble + 4.25) * t_sym
    
    # Payload time
    t_payload = pl_sym * t_sym
    
    return t_preamble + t_payload

def main():
    # --- Parameters ---
    payload = 100  # bytes
    sf_options = range(7, 13)
    tx_power_options = range(2, 15, 2)  # 2 to 14 dBm in 2dB steps
    
    # SNR demodulation thresholds for each SF (dB)
    snr_thresholds = {
        7: -7.5,
        8: -10.0,
        9: -12.5,
        10: -15.0,
        11: -17.5,
        12: -20.0
    }
    
    # Fading margin for Rician K=3dB at 1% PER (99% reliability)
    fading_margin = 4.8  # dB
    
    # Assumed total link loss for an urban environment (Path Loss + Noise Figure - Gains)
    assumed_link_loss = 135.0  # dB
    
    # Thermal noise floor adjustment term (-174dBm/Hz + 10*log10(BW))
    # For BW=125kHz, this is approx. -120 dBm. So P_rx = P_tx - LinkLoss + 120 + SNR
    noise_floor_adjustment = 120.0 # dB

    results = []
    
    for sf in sf_options:
        # 1. Calculate Time on Air
        toa = calculate_toa(payload, sf)
        
        # 2. Determine Required SNR
        snr_th = snr_thresholds[sf]
        snr_avg_req = snr_th + fading_margin
        
        # 3. Calculate theoretically required transmit power for the assumed link loss
        # P_req = LinkLoss - NoiseFloorAdjustment + SNR_avg_req
        p_req_dbm = assumed_link_loss - noise_floor_adjustment + snr_avg_req
        
        # 4. Find the lowest available transmit power that meets the requirement
        chosen_tx_power_dbm = -1
        if p_req_dbm > max(tx_power_options):
            # This SF cannot meet the link budget even at max power
            continue 
            
        for p in tx_power_options:
            if p >= p_req_dbm:
                chosen_tx_power_dbm = p
                break
        
        # 5. Calculate energy consumption for this configuration
        # Energy is proportional to Power (mW) * Time (s)
        power_mw = 10**(chosen_tx_power_dbm / 10)
        energy = power_mw * toa
        
        results.append({
            "sf": sf,
            "p_req_dbm": p_req_dbm,
            "chosen_tx_power_dbm": chosen_tx_power_dbm,
            "energy": energy,
            "snr_avg_req": snr_avg_req
        })

    # 6. Find the optimal configuration with the minimum energy
    if not results:
        print("No configuration could meet the required link budget.")
        return
        
    optimal_config = min(results, key=lambda x: x['energy'])
    
    # --- Final Output ---
    print("Based on an assumed urban link loss of {} dB:".format(assumed_link_loss))
    print("\nThe optimal Spreading Factor is: SF{}".format(optimal_config['sf']))
    print("The optimal Transmission Power is: {} dBm".format(optimal_config['chosen_tx_power_dbm']))
    
    print("\nThis is derived from the required transmission power calculation:")
    
    print("Required Tx Power (dBm) = Link Loss (dB) - Noise Floor Adjustment (dB) + Required Average SNR (dB)")
    
    print(f"Required Tx Power (dBm) = {assumed_link_loss} - {noise_floor_adjustment} + ({optimal_config['snr_avg_req']:.1f}) = {optimal_config['p_req_dbm']:.1f} dBm")

if __name__ == '__main__':
    main()
