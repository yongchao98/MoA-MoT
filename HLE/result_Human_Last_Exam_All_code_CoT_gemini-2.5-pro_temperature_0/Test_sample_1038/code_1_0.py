import math

def calculate_resistance_improvement():
    """
    This function models the effect of body biasing on a subthreshold transistor
    as described in option C. It calculates the ratio of resistance in the
    biased (operate) state versus the unbiased (reset) state.
    """
    # --- Parameters ---
    # Assume a gate-source voltage that is slightly below the threshold
    v_gs = 0.4  # Volts
    
    # Nominal threshold voltage (Vt) of the transistor
    vt_nominal = 0.45  # Volts
    
    # Subthreshold slope factor (n), typically 1.3-1.7
    n = 1.5
    
    # Thermal voltage (kT/q) at room temperature
    U_T = 0.026  # Volts
    
    # The body bias in 'operate' mode increases the effective threshold voltage.
    # Let's assume it increases Vt by 0.2V.
    vt_increase_from_bias = 0.2
    vt_biased = vt_nominal + vt_increase_from_bias

    # The subthreshold current is proportional to exp((Vgs - Vt) / (n * U_T))
    # We can ignore the constant pre-factor I_0 as we are interested in the ratio.
    
    # --- Calculations ---
    
    # 1. Calculate exponent for the nominal case (lower resistance for reset)
    exponent_nominal = (v_gs - vt_nominal) / (n * U_T)
    # The relative current is proportional to math.exp(exponent_nominal)
    relative_current_nominal = math.exp(exponent_nominal)

    # 2. Calculate exponent for the body-biased case (higher resistance for operate)
    exponent_biased = (v_gs - vt_biased) / (n * U_T)
    # The relative current is proportional to math.exp(exponent_biased)
    relative_current_biased = math.exp(exponent_biased)
    
    # 3. Resistance is inversely proportional to current (R ~ V/I).
    # The ratio of resistances is the inverse of the ratio of currents.
    resistance_improvement_factor = relative_current_nominal / relative_current_biased

    # --- Output ---
    print("--- Analysis of Body Biasing (Option C) ---")
    print(f"Nominal Threshold Voltage (Vt): {vt_nominal} V")
    print(f"Body-Biased Threshold Voltage (Vt'): {vt_biased:.2f} V")
    print(f"Assumed Gate-Source Voltage (Vgs): {v_gs} V")
    print("\nSubthreshold current is proportional to exp((Vgs - Vt) / (n * U_T))")
    print(f"Exponent for nominal Vt: ({v_gs} - {vt_nominal}) / ({n} * {U_T}) = {exponent_nominal:.2f}")
    print(f"Exponent for biased Vt': ({v_gs} - {vt_biased}) / ({n} * {U_T}) = {exponent_biased:.2f}")
    
    print(f"\nRelative subthreshold current (nominal): exp({exponent_nominal:.2f}) = {relative_current_nominal:.4e}")
    print(f"Relative subthreshold current (biased): exp({exponent_biased:.2f}) = {relative_current_biased:.4e}")
    
    print("\n--- Result ---")
    print("Since resistance is inversely proportional to current, the resistance in the 'operate' mode is significantly higher.")
    print(f"Resistance Improvement Factor: {resistance_improvement_factor:,.2f} times higher")
    print("This demonstrates that body biasing is highly effective at creating the high resistance needed for the 'operate' phase.")

calculate_resistance_improvement()
<<<C>>>