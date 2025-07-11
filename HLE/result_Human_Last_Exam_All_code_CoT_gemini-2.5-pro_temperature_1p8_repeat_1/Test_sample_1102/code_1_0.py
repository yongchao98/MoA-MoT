def find_welding_parameters():
    """
    This function analyzes the welding parameters to determine the appropriate current and voltage.

    Key considerations:
    1. Arc Gap and Voltage: A long arc gap of 6 mm requires a correspondingly high voltage to sustain the arc.
       Voltages below ~12V are generally not viable for such a long arc in TIG welding.
    2. Material and Heat Input: Inconel 718 blade tips are thin and sensitive to heat.
       A slow travel speed (0.5 mm/s) increases heat input per unit length.
       To prevent overheating, distortion, or burn-through, a low current (amperage) is necessary.

    Evaluation of choices leads to the selection of a low current and a high voltage.
    The chosen parameters provide a stable arc (due to high voltage) with controlled, low heat input (due to low current),
    which is ideal for the repair of a thin Inconel 718 turbine blade.
    """
    # Based on the analysis, the most suitable parameters are a low current and a high voltage.
    current_amps = 25
    voltage_volts = 17.5

    # The chosen answer corresponds to these values.
    # The output format is the final equation with each number printed.
    print(f"The recommended welding specification is: Current = {current_amps} A and Voltage = {voltage_volts} V")

find_welding_parameters()