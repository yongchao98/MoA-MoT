import math
from scipy.special import marcumq

def solve_lorawan_optimization():
    """
    Calculates the optimal LoRaWAN Spreading Factor and Transmit Power
    to minimize energy while maintaining a target Packet Error Rate (PER).
    """
    # --- Given Parameters ---
    payload_bytes = 100
    bandwidth_hz = 125000
    coding_rate_numerator = 4
    coding_rate_denominator = 5
    k_factor_db = 3
    target_per = 0.01
    transmit_powers_dbm = list(range(2, 15, 2))
    spreading_factors = list(range(7, 13))

    # --- Model Assumptions & Constants ---
    # Path Loss Calculation (Okumura-Hata for Urban, 1km distance)
    # This is a reasonable assumption for a fixed urban link.
    # The absolute value isn't as critical as its consistent application.
    path_loss_db = 126

    # Receiver Noise Floor
    noise_figure_db = 6 # Typical gateway receiver noise figure
    noise_floor_dbm = -174 + 10 * math.log10(bandwidth_hz) + noise_figure_db

    # LoRaWAN Receiver Sensitivity (minimum required SNR for demodulation)
    snr_req_db = {
        7: -7.5, 8: -10.0, 9: -12.5,
        10: -15.0, 11: -17.5, 12: -20.0
    }

    # Rician K-factor in linear scale
    k_factor_linear = 10**(k_factor_db / 10)

    # LoRaWAN constants for Time on Air calculation
    preamble_syms = 8
    header_enabled = 0  # 0 for explicit header
    crc_enabled = 1     # 1 for payload CRC on

    # --- Optimization ---
    optimal_settings = None
    min_energy_joules = float('inf')

    print("--- LoRaWAN ADR Optimization Analysis ---")
    print(f"Target PER: < {target_per * 100}%\n")

    # Iterate through each Spreading Factor
    for sf in spreading_factors:
        # 1. Calculate Time on Air (ToA) for this SF
        symbol_time_s = (2**sf) / bandwidth_hz
        low_dr_optimize = 1 if symbol_time_s > 0.016 else 0
        
        # Numerator of the payload symbols calculation
        payload_syms_num = (8 * payload_bytes - 4 * sf + 28 + 16 * crc_enabled - 20 * header_enabled)
        # Denominator of the payload symbols calculation
        payload_syms_den = 4 * (sf - 2 * low_dr_optimize)
        
        # Coded payload symbols
        payload_syms = math.ceil(max(0, payload_syms_num) / payload_syms_den) * (coding_rate_denominator)

        # Total symbols in the packet
        total_syms = preamble_syms + 4.25 + 8 + payload_syms
        time_on_air_s = total_syms * symbol_time_s

        # 2. Find minimum power for this SF that meets PER target
        for tx_power_dbm in transmit_powers_dbm:
            # Calculate average SNR at the receiver
            avg_snr_db = tx_power_dbm - path_loss_db - noise_floor_dbm
            avg_snr_linear = 10**(avg_snr_db / 10)

            # Get required SNR for this SF
            req_snr_linear = 10**(snr_req_db[sf] / 10)

            # Calculate PER using Rician fading model (outage probability)
            # PER = 1 - Q1(a, b), where Q1 is the Marcum Q-function of order 1
            a = math.sqrt(2 * k_factor_linear)
            b = math.sqrt(2 * (k_factor_linear + 1) * (req_snr_linear / avg_snr_linear))
            
            # scipy.special.marcumq(a, b, 1) is the Marcum Q-function Q_1(a,b)
            calculated_per = 1.0 - marcumq(a, b, 1)

            # Check if the PER requirement is met
            if calculated_per <= target_per:
                # This is a valid setting. Calculate energy.
                tx_power_watts = 10**((tx_power_dbm - 30) / 10)
                energy_joules = tx_power_watts * time_on_air_s

                # Check if this is the best solution so far
                if energy_joules < min_energy_joules:
                    min_energy_joules = energy_joules
                    optimal_settings = {
                        "SF": sf,
                        "TxPower": tx_power_dbm,
                        "Energy": energy_joules,
                        "ToA": time_on_air_s,
                        "PER": calculated_per
                    }
                # Since we iterate power upwards, this is the lowest power for this SF.
                # No need to check higher powers for this SF.
                break

    # --- Print Final Result ---
    if optimal_settings:
        print("Optimal settings found:")
        print(f"  Spreading Factor (SF): {optimal_settings['SF']}")
        print(f"  Transmission Power (TxP): {optimal_settings['TxPower']} dBm")
        print("\nDerived Metrics for this setting:")
        print(f"  Time on Air (ToA): {optimal_settings['ToA']:.3f} seconds")
        print(f"  Energy per Transmission: {optimal_settings['Energy'] * 1000:.3f} mJ")
        print(f"  Calculated Packet Error Rate (PER): {optimal_settings['PER']:.4f} ({optimal_settings['PER']*100:.2f}%)")
        
        # Final answer format
        return f"<<<SF{optimal_settings['SF']}, {optimal_settings['TxPower']} dBm>>>"
    else:
        print("No combination of SF and TxPower could meet the PER requirement with the given model.")
        return "<<<No solution found>>>"

# Execute the function and print the result
final_answer = solve_lorawan_optimization()
print(final_answer)