import math

def solve_lorawan_optimization():
    """
    Calculates the optimal LoRaWAN Spreading Factor and Transmit Power
    to minimize energy consumption for a given scenario.
    """

    # 1. Define physical and LoRaWAN parameters
    PL = 100  # Payload in bytes
    BW = 125000  # Bandwidth in Hz
    CR_CODE = 1  # Coding Rate 4/5 corresponds to CR_code=1
    N_PREAMBLE = 8 # Standard number of preamble symbols

    # Available hardware settings
    TX_POWERS_DBM = list(range(2, 16, 2))  # 2, 4, 6, 8, 10, 12, 14 dBm
    SPREADING_FACTORS = list(range(7, 13)) # SF7 to SF12

    # Receiver sensitivity (minimum required SNR for demodulation) for each SF
    SNR_REQUIREMENTS = {
        7: -7.5,
        8: -10.0,
        9: -12.5,
        10: -15.0,
        11: -17.5,
        12: -20.0,
    }

    # 2. Define channel and environmental assumptions
    # Based on a Rician fading channel with K=3dB and 99% reliability (1% PER)
    FADING_MARGIN_DB = 5.5
    # A representative path loss for an urban LoRaWAN link
    PATH_LOSS_DB = 138.0
    # A typical noise figure for a LoRaWAN gateway receiver
    NOISE_FIGURE_DB = 6.0

    # 3. Calculate constant values
    # Noise Power (Pn) = NoiseFigure + 10*log10(BW) + Thermal Noise Floor
    # Thermal Noise Floor at room temperature is -174 dBm/Hz
    NOISE_POWER_DBM = NOISE_FIGURE_DB + 10 * math.log10(BW) - 174

    print("--- System Assumptions ---")
    print(f"Path Loss: {PATH_LOSS_DB} dB")
    print(f"Fading Margin (Rician K=3dB, PER<1%): {FADING_MARGIN_DB} dB")
    print(f"Receiver Noise Power: {NOISE_POWER_DBM:.2f} dBm\n")


    def calculate_toa_s(sf):
        """Calculates Time on Air (ToA) in seconds for a given SF."""
        # Low Data Rate Optimization (LDRO) is enabled for SF > 10 on 125kHz
        de = 1 if sf > 10 else 0
        
        # Symbol time in seconds
        t_sym_s = (2**sf) / BW

        # Preamble duration
        t_preamble_s = (N_PREAMBLE + 4.25) * t_sym_s

        # Payload calculation
        # Using standard formula for explicit header (H=0)
        payload_numerator = 8 * PL - 4 * sf + 28 + 16
        payload_denominator = 4 * (sf - 2 * de)
        
        # Number of symbols in the payload
        n_payload_sym = 8 + max(0, math.ceil(payload_numerator / payload_denominator) * (CR_CODE + 4))

        # Payload duration
        t_payload_s = n_payload_sym * t_sym_s
        
        return t_preamble_s + t_payload_s

    # 4. Optimization Loop
    min_energy_joules = float('inf')
    optimal_params = None

    print("--- Analyzing Viable (SF, Power) Combinations ---")
    for sf in SPREADING_FACTORS:
        # Calculate the total required signal power at the transmitter
        snr_req = SNR_REQUIREMENTS[sf]
        required_txp_dbm = snr_req + PATH_LOSS_DB + NOISE_POWER_DBM + FADING_MARGIN_DB

        # Find the minimum available power that meets the requirement
        min_viable_txp = None
        for txp in TX_POWERS_DBM:
            if txp >= required_txp_dbm:
                min_viable_txp = txp
                break
        
        # If no power level is sufficient, this SF is not viable
        if min_viable_txp is None:
            print(f"SF{sf}: Not viable (Required TxP > {TX_POWERS_DBM[-1]} dBm)")
            continue

        # Calculate energy for this viable combination
        toa_s = calculate_toa_s(sf)
        power_watts = 10**((min_viable_txp - 30) / 10)
        energy_joules = power_watts * toa_s

        print(f"SF{sf}: Viable with {min_viable_txp} dBm. ToA={toa_s:.3f} s, Energy={energy_joules*1000:.2f} mJ")

        # Check if this is the new most energy-efficient setting
        if energy_joules < min_energy_joules:
            min_energy_joules = energy_joules
            optimal_params = {
                "sf": sf,
                "txp_dbm": min_viable_txp,
                "power_watts": power_watts,
                "toa_s": toa_s,
                "energy_j": energy_joules
            }

    # 5. Output the final result
    print("\n--- Optimal Configuration ---")
    if optimal_params:
        sf = optimal_params["sf"]
        txp = optimal_params["txp_dbm"]
        p_w = optimal_params["power_watts"]
        toa = optimal_params["toa_s"]
        energy = optimal_params["energy_j"]
        
        print(f"The optimal setting to minimize energy while keeping PER < 1% is:")
        print(f"Spreading Factor: SF{sf}")
        print(f"Transmission Power: {txp} dBm")
        print("\nFinal Energy Calculation:")
        # The prompt requires printing the numbers in the final equation.
        print(f"Energy (J) = Transmit Power (W) * Time on Air (s)")
        print(f"Energy (J) = {p_w:.4f} * {toa:.4f} = {energy:.4f} J ({energy*1000:.2f} mJ)")
        
        # Final answer format as requested
        print(f"\n<<<SF={sf}, TxP={txp} dBm>>>")

    else:
        print("No viable combination of SF and TxP could be found for the given path loss and PER requirement.")
        print("\n<<<No solution found>>>")

# Execute the solver
solve_lorawan_optimization()