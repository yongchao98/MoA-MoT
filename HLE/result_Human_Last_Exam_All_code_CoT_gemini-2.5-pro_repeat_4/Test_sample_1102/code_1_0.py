# Based on the analysis of welding physics for the given parameters.

# 1. Voltage requirement for a 6 mm TIG arc gap.
# A stable TIG arc requires a voltage of approximately 10-12V plus 1V per mm of arc length.
# For a 6 mm arc gap, the required voltage is ~11V + (1V/mm * 6mm) = 17V.
# The closest and most plausible voltage from the choices is 17.5 V.

# 2. Heat input control for Inconel 718 at a slow travel speed.
# The material is sensitive to high heat input, and the travel speed is very slow (0.5 mm/s).
# To avoid overheating and cracking, a low welding current is necessary.
# 25 A is a suitable low current for a controlled root pass build-up under these conditions.
# Heat Input (HI) = (Voltage * Current) / Travel Speed
# HI = (17.5 V * 25 A) / 0.5 mm/s = 875 J/mm, which is a controlled and moderate value.

current_amps = 25
voltage_volts = 17.5

print("Based on the physics of TIG welding for the specified conditions, the recommended parameters are:")
print(f"Current: {current_amps} A")
print(f"Voltage: {voltage_volts} V")