import math

def analyze_body_bias_effect():
    """
    Analyzes the effect of body biasing on subthreshold leakage current as
    described in option C.

    The subthreshold current (I_sub) is proportional to exp((Vgs - Vt) / (n * V_thermal)),
    where Vt is the threshold voltage. By increasing Vt, we can exponentially
    decrease the current for a given Vgs.
    """

    # --- Given Parameters ---
    vt_native = 0.45  # volts, native threshold voltage
    # Option C suggests a body bias that effectively raises Vt.
    # The increase in Vt due to body effect can be significant.
    # Let's assume the body bias effectively increases Vt by 0.25V.
    vt_increase = 0.25 # volts, an example effective increase in Vt from body bias
    vt_body_biased = vt_native + vt_increase

    # Subthreshold slope factor (n), let's use the average of the given range [1.4, 1.6]
    n = 1.5

    # Thermal voltage at room temperature (300K)
    k = 1.380649e-23  # Boltzmann constant in J/K
    T = 300           # Temperature in Kelvin
    q = 1.60217663e-19 # Elementary charge in Coulombs
    V_thermal = k * T / q  # Thermal voltage in volts, approximately 26mV

    # --- Calculation ---
    # The term in the denominator of the exponent
    m = n * V_thermal

    # The ratio of subthreshold current with body bias to without body bias is:
    # Ratio = exp((Vgs - vt_body_biased) / m) / exp((Vgs - vt_native) / m)
    # This simplifies to exp(- (vt_body_biased - vt_native) / m), which is exp(-vt_increase / m)
    current_reduction_factor = math.exp(-vt_increase / m)

    # --- Output Results ---
    print("--- Analysis of Design Strategy C: Body Biasing ---")
    print(f"Native Threshold Voltage (vt_native): {vt_native:.2f} V")
    print(f"Assumed Vt Increase from Body Bias: {vt_increase:.2f} V")
    print(f"New Effective Threshold Voltage (vt_body_biased): {vt_body_biased:.2f} V")
    print("\nCalculating the subthreshold current reduction:")
    print(f"Subthreshold Slope Factor (n): {n}")
    print(f"Thermal Voltage (V_thermal): {V_thermal:.4f} V")
    print(f"Denominator term (n * V_thermal): {m:.4f} V")
    print("\nEquation for reduction: exp(-vt_increase / (n * V_thermal))")
    print(f"Plugging in numbers: exp(-{vt_increase:.2f} / ({n} * {V_thermal:.4f})) = exp(-{vt_increase:.2f} / {m:.4f})")
    print(f"\nResult:")
    print(f"The subthreshold leakage current is reduced by a factor of approximately {1/current_reduction_factor:.0f}.")
    print("\nThis massive reduction in leakage current is why option C is the most effective strategy.")
    print("It allows for a large, stable gate capacitor to be used (solving the leakage and charge injection problem) while still achieving very high resistance during operation.")

analyze_body_bias_effect()