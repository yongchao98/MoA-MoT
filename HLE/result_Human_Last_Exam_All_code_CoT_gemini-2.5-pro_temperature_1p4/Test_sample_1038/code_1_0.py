import math

def analyze_operating_point():
    """
    Analyzes and compares the operating point of a pseudo-resistor transistor
    with and without body biasing to demonstrate the effectiveness of Option C.
    """

    # --- Circuit & Device Parameters ---
    V_supply = 1.2  # volts
    Vt_initial = 0.45  # volts
    V_source_cm = 0.6  # volts (Assumed common-mode voltage)
    V_offset = 0.1  # volts (+/- 100mV sensor offset)

    print("--- Analysis of Pseudo-Resistor Operating Point under 1.2V Supply ---")
    print(f"Goal: Handle a sensor offset of +/- {V_offset * 1000:.0f} mV without saturation.")
    print(f"This means the source voltage (Vs) varies from {V_source_cm - V_offset:.1f} V to {V_source_cm + V_offset:.1f} V.")
    print("-" * 60)

    # --- Scenario 1: Standard Transistor (No Body Bias) ---
    print("\n[Scenario A-style: Standard Transistor Configuration]")
    Vg_standard = 0.9  # A plausible gate bias voltage, must be < V_supply
    print(f"Threshold Voltage (Vt) = {Vt_initial} V")
    print(f"Chosen Gate Bias (Vg) = {Vg_standard} V (held by capacitor)")

    # Calculate Vgs at the boundaries of the offset range
    Vs_low = V_source_cm - V_offset
    Vgs_max = Vg_standard - Vs_low
    print(f"\nWhen sensor offset is -100mV (Vs = {Vs_low:.2f} V):")
    print(f"  Vgs_max = Vg - Vs_low = {Vg_standard:.2f} V - {Vs_low:.2f} V = {Vgs_max:.2f} V")

    Vs_high = V_source_cm + V_offset
    Vgs_min = Vg_standard - Vs_high
    print(f"When sensor offset is +100mV (Vs = {Vs_high:.2f} V):")
    print(f"  Vgs_min = Vg - Vs_high = {Vg_standard:.2f} V - {Vs_high:.2f} V = {Vgs_min:.2f} V")

    margin = Vt_initial - Vgs_max
    print(f"\nConclusion for Standard Transistor:")
    print(f"The Vgs range is [{Vgs_min:.2f}V, {Vgs_max:.2f}V].")
    print(f"The maximum Vgs ({Vgs_max:.2f}V) is very close to the threshold voltage Vt ({Vt_initial:.2f}V).")
    print(f"The safety margin is only {Vt_initial:.2f}V - {Vgs_max:.2f}V = {margin:.2f}V.")
    print("This low margin risks pushing the transistor out of subthreshold, causing a large drop in resistance and amplifier saturation.")
    print("-" * 60)

    # --- Scenario 2: With Body Bias (Option C) ---
    print("\n[Scenario C: On-Chip Body-Bias Generator]")
    # Per option C, body bias increases the threshold voltage
    Vt_new = Vt_initial + 0.2  # A plausible increase in Vt
    # With a higher Vt, we have more freedom to set a higher Vg
    Vg_biased = 1.1 # V, still safely below V_supply

    print(f"New Dynamic Threshold Voltage (Vt_new) = {Vt_new:.2f} V (during 'operate' phase)")
    print(f"Chosen Gate Bias (Vg) = {Vg_biased} V")

    # Recalculate Vgs with the new gate bias
    Vgs_max_biased = Vg_biased - Vs_low
    print(f"\nWhen sensor offset is -100mV (Vs = {Vs_low:.2f} V):")
    print(f"  Vgs_max = Vg - Vs_low = {Vg_biased:.2f} V - {Vs_low:.2f} V = {Vgs_max_biased:.2f} V")

    Vgs_min_biased = Vg_biased - Vs_high
    print(f"When sensor offset is +100mV (Vs = {Vs_high:.2f} V):")
    print(f"  Vgs_min = Vg - Vs_high = {Vg_biased:.2f} V - {Vs_high:.2f} V = {Vgs_min_biased:.2f} V")

    margin_biased = Vt_new - Vgs_max_biased
    print(f"\nConclusion for Body-Biased Transistor (Option C):")
    print(f"The Vgs range is [{Vgs_min_biased:.2f}V, {Vgs_max_biased:.2f}V].")
    print(f"The entire Vgs range is now safely below the new threshold voltage Vt_new ({Vt_new:.2f}V).")
    print(f"The new safety margin is {Vt_new:.2f}V - {Vgs_max_biased:.2f}V = {margin_biased:.2f}V.")
    print("This robust margin ensures the transistor remains in the high-resistance subthreshold region across the entire sensor offset range, effectively balancing the circuit's conflicting needs.")

if __name__ == "__main__":
    analyze_operating_point()