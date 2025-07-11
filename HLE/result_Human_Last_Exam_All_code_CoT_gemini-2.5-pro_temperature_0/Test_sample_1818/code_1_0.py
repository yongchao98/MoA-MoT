import numpy as np
from scipy.special import marcumq
from scipy.optimize import brentq
import math

def solve_lorawan_optimization():
    """
    Calculates the optimal Spreading Factor and Transmit Power for a LoRaWAN device
    to minimize energy consumption while meeting a PER requirement.
    """
    # --- 1. Constants and Parameters ---
    # LoRaWAN Parameters
    PAYLOAD = 100  # bytes
    BW = 125e3  # Hz
    CR_CODE = 1 # Coding Rate is 4/5 -> 4/(4+1)
    PREAMBLE_LEN = 8
    HEADER_ENABLED = True # Explicit header
    CRC_ENABLED = True
    LOW_DR_OPTIMIZE = True # As per LoRaWAN spec for SF11, SF12

    # Channel and Link Parameters
    K_FACTOR_DB = 3
    K_FACTOR_LIN = 10**(K_FACTOR_DB / 10)
    TARGET_PER = 0.01
    # Assumption: A typical path loss for an urban environment. The optimal TxP depends on this value.
    PATH_LOSS_DB = 135
    # Typical gateway noise figure
    NOISE_FIGURE_DB = 6
    THERMAL_NOISE_DBM = -174 # dBm/Hz
    NOISE_FLOOR_DBM = THERMAL_NOISE_DBM + 10 * np.log10(BW) + NOISE_FIGURE_DB

    # Available Device Settings
    TX_POWERS_DBM = list(range(2, 15, 2))
    SPREADING_FACTORS = list(range(7, 13))

    # SNR thresholds for demodulation (from Semtech datasheets)
    SNR_THRESHOLDS_DB = {
        7: -7.5, 8: -10, 9: -12.5, 10: -15, 11: -17.5, 12: -20
    }

    # --- 2. Helper Functions ---

    def calculate_required_snr(snr_threshold_db, k_factor_lin, target_per):
        """Calculates the required average SNR on a Rician channel."""
        snr_th_lin = 10**(snr_threshold_db / 10)
        a = np.sqrt(2 * k_factor_lin)
        
        def equation_to_solve(snr_avg_lin):
            if snr_avg_lin <= 0: return 1
            b = np.sqrt(2 * (k_factor_lin + 1) * snr_th_lin / snr_avg_lin)
            # PER = 1 - MarcumQ(a, b)
            return (1 - marcumq(a, b)) - target_per

        try:
            # The average SNR must be higher than the threshold. We search in a reasonable range.
            required_snr_lin = brentq(equation_to_solve, snr_th_lin, 1000 * snr_th_lin)
            return 10 * np.log10(required_snr_lin)
        except ValueError:
            return float('inf')

    def calculate_toa(sf, payload, bw, cr_code):
        """Calculates the Time on Air (ToA) for a LoRa packet."""
        de = 1 if sf in [11, 12] else 0
        t_sym = (2**sf) / bw
        t_preamble = (PREAMBLE_LEN + 4.25) * t_sym
        h = 0 # Header enabled
        crc = 16 # CRC enabled
        
        payload_symb_nb_num = 8 * payload - 4 * sf + 28 + crc - 20 * h
        payload_symb_nb_den = 4 * (sf - 2 * de)
        
        payload_symb_nb = 8 + max(0, math.ceil(payload_symb_nb_num / payload_symb_nb_den) * (cr_code + 4))
        t_payload = payload_symb_nb * t_sym
        return t_preamble + t_payload

    # --- 3. Main Calculation Loop ---
    results = []
    min_energy = float('inf')
    optimal_config = None

    print("Calculating optimal LoRaWAN parameters...")
    print(f"Assumptions: Path Loss = {PATH_LOSS_DB} dB, Gateway Noise Figure = {NOISE_FIGURE_DB} dB")
    print("-" * 80)
    print(f"{'SF':<5} {'SNR_req (dB)':<15} {'ToA (ms)':<12} {'TxP_req (dBm)':<15} {'Min TxP (dBm)':<15} {'Energy (mJ)':<15}")
    print("-" * 80)

    for sf in SPREADING_FACTORS:
        required_snr_db = calculate_required_snr(SNR_THRESHOLDS_DB[sf], K_FACTOR_LIN, TARGET_PER)
        toa_s = calculate_toa(sf, PAYLOAD, BW, CR_CODE)
        
        # Required Tx Power = Required SNR at receiver + Path Loss + Noise Floor
        required_tx_power_dbm = required_snr_db + PATH_LOSS_DB + NOISE_FLOOR_DBM
        
        min_tx_power_dbm = float('inf')
        for p in TX_POWERS_DBM:
            if p >= required_tx_power_dbm:
                min_tx_power_dbm = p
                break
            
        energy_mj = float('inf')
        if min_tx_power_dbm != float('inf'):
            tx_power_w = 10**((min_tx_power_dbm - 30) / 10)
            energy_j = tx_power_w * toa_s
            energy_mj = energy_j * 1000
            
            if energy_mj < min_energy:
                min_energy = energy_mj
                optimal_config = {
                    'sf': sf, 'tx_power': min_tx_power_dbm, 'energy': energy_mj,
                    'toa': toa_s, 'req_snr': required_snr_db
                }

        print(f"SF{sf:<2} {required_snr_db:<15.2f} {toa_s*1000:<12.2f} {required_tx_power_dbm:<15.2f} {str(min_tx_power_dbm):<15} {energy_mj:<15.2f}")

    # --- 4. Final Answer ---
    print("-" * 80)
    if optimal_config:
        print("\nOptimal Configuration Found:")
        print(f"  Spreading Factor: SF{optimal_config['sf']}")
        print(f"  Transmit Power:   {optimal_config['tx_power']} dBm")
        print(f"  Lowest Energy:    {optimal_config['energy']:.2f} mJ")
        
        sf_opt = optimal_config['sf']
        txp_opt = optimal_config['tx_power']
        toa_opt = optimal_config['toa']
        
        print("\nFinal Energy Calculation for the Optimal Configuration:")
        print(f"Energy (mJ) = 10^(({txp_opt} - 30) / 10) * {toa_opt:.4f} * 1000 = {optimal_config['energy']:.2f}")
    else:
        print("\nNo viable configuration found for the assumed path loss.")
        print("The required transmit power exceeds the maximum available (14 dBm) for all SFs.")

solve_lorawan_optimization()