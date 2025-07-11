import math
import numpy as np
from scipy.special import marcumq
from scipy.optimize import root_scalar

def solve_lorawan_optimization():
    """
    Calculates the optimal LoRaWAN Spreading Factor and Transmit Power
    to minimize energy consumption for a given set of constraints.
    """

    # --- Step 1: Define Given Parameters ---
    PAYLOAD_BYTES = 100
    K_FACTOR_DB = 3
    PER_TARGET = 0.01
    BANDWIDTH_HZ = 125000
    CODING_RATE_CODE = 1  # For CR 4/5
    N_PREAMBLE = 8
    HEADER_ENABLED = True
    CRC_ENABLED = True
    AVAILABLE_TP_DBM = list(range(2, 15, 2))
    AVAILABLE_SF = list(range(7, 13))

    # SNR demodulation thresholds for each Spreading Factor (from datasheets)
    SNR_THRESHOLDS_DB = {
        7: -7.5, 8: -10, 9: -12.5, 10: -15, 11: -17.5, 12: -20
    }

    # --- Step 2: Assume Path Loss and Calculate Noise Floor ---
    # We assume a representative path loss for an urban environment.
    # This value is critical for determining the required transmit power.
    PATH_LOSS_DB = 136
    
    # Standard thermal noise at room temperature for the given bandwidth
    NOISE_FLOOR_DBM = -174 + 10 * math.log10(BANDWIDTH_HZ)
    # A typical noise figure for a gateway receiver
    GATEWAY_NOISE_FIGURE_DB = 6
    TOTAL_NOISE_DBM = NOISE_FLOOR_DBM + GATEWAY_NOISE_FIGURE_DB

    print(f"Analyzing for an assumed Path Loss of {PATH_LOSS_DB} dB...\n")

    # --- Helper Functions ---
    def find_marcum_b(a, q_target):
        """Finds b such that MarcumQ(a, b, 1) = q_target using a root finder."""
        def func_to_solve(b):
            return marcumq(a, b, 1) - q_target
        # b must be non-negative. Bracket is chosen to be wide enough.
        sol = root_scalar(func_to_solve, bracket=[0, 50], method='brentq')
        return sol.root

    def calculate_toa_s(sf):
        """Calculates the Time on Air (ToA) in seconds."""
        de = 1 if sf >= 11 else 0
        h = 0 if HEADER_ENABLED else 1
        crc_val = 16 if CRC_ENABLED else 0
        
        t_sym = (2**sf) / BANDWIDTH_HZ
        t_preamble = (n_preamble + 4.25) * t_sym
        
        numerator = 8 * PAYLOAD_BYTES - 4 * sf + 28 + crc_val - 20 * h
        denominator = 4 * (sf - 2 * de)
        
        payload_symb_nb = 8 + max(math.ceil(numerator / denominator) * (CODING_RATE_CODE + 4), 0)
        t_payload = payload_symb_nb * t_sym
        return t_preamble + t_payload

    # --- Step 3: Determine Required SNR for PER <= 1% in Rician Channel ---
    K_LINEAR = 10**(K_FACTOR_DB / 10)
    a_marcum = math.sqrt(2 * K_LINEAR)
    # We need PER <= 0.01, so we solve for the boundary case PER = 0.01
    # PER = 1 - Q1(a, b), so Q1(a, b) = 1 - PER = 0.99
    b_target = find_marcum_b(a_marcum, 1 - PER_TARGET)

    required_avg_snr_db = {}
    for sf in AVAILABLE_SF:
        snr_req_linear = 10**(SNR_THRESHOLDS_DB[sf] / 10)
        # From b = sqrt(2*(K+1)*snr_req/snr_avg), we solve for snr_avg
        snr_avg_min_linear = (2 * (K_LINEAR + 1) * snr_req_linear) / (b_target**2)
        required_avg_snr_db[sf] = 10 * math.log10(snr_avg_min_linear)

    # --- Step 4 & 5: Iterate, Calculate Energy, and Find Optimum ---
    min_energy = float('inf')
    optimal_sf = None
    optimal_tp = None
    
    print("--- Energy Consumption Analysis ---")
    
    for sf in AVAILABLE_SF:
        # Calculate the TP required to overcome path loss and achieve the target SNR
        required_tp_dbm = required_avg_snr_db[sf] - TOTAL_NOISE_DBM + PATH_LOSS_DBM
        
        # Find the smallest available TP that meets the requirement
        chosen_tp_dbm = None
        for tp in AVAILABLE_TP_DBM:
            if tp >= required_tp_dbm:
                chosen_tp_dbm = tp
                break
        
        # If a valid TP is found, calculate energy and compare
        if chosen_tp_dbm is not None:
            toa_s = calculate_toa_s(sf)
            tp_mw = 10**(chosen_tp_dbm / 10)
            
            # Energy is proportional to Power * Time on Air.
            # We use units of mW * ms for clarity.
            energy = tp_mw * (toa_s * 1000)
            
            print(f"SF{sf}:")
            print(f"  - Required Avg. SNR: {required_avg_snr_db[sf]:.2f} dB")
            print(f"  - Required TP: {required_tp_dbm:.2f} dBm -> Chosen TP: {chosen_tp_dbm} dBm")
            print(f"  - Time on Air: {toa_s*1000:.1f} ms")
            print(f"  - Energy Calculation: {tp_mw:.1f} mW * {toa_s*1000:.1f} ms = {energy:.0f} (relative units)")
            
            if energy < min_energy:
                min_energy = energy
                optimal_sf = sf
                optimal_tp = chosen_tp_dbm
        else:
            print(f"SF{sf}:")
            print(f"  - Required TP ({required_tp_dbm:.2f} dBm) exceeds maximum available power (14 dBm). Not a viable option.")
        print("-" * 35)

    if optimal_sf is not None:
        print("\n--- Optimal Configuration ---")
        print(f"The optimal Spreading Factor is SF{optimal_sf}.")
        print(f"The optimal Transmission Power is {optimal_tp} dBm.")
    else:
        print("\nNo viable configuration found for the assumed path loss.")

solve_lorawan_optimization()
<<<SF10, 10 dBm>>>