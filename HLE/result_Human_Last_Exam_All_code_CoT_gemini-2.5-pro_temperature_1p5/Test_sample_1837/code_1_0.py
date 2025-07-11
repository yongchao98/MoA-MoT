import math

# Step 1: Define constants and the user's measured values
h = 6.62607015e-34  # Planck's constant (J·s)
c = 299792458      # Speed of light (m/s)
k = 1.380649e-23   # Boltzmann constant (J/K)

L_nm = 400         # Wavelength in nanometers
T_measured = 9000  # Surface temperature in Kelvin
B_measured = 1.2e15 # Measured spectral radiance in W·sr⁻¹·m⁻²·m⁻¹

# Convert wavelength to meters for calculation
L = L_nm * 1e-9

# Step 2: Explain the physical law and the plan
print("The relationship between a black body's spectral radiance (B), temperature (T), and wavelength (L) is given by Planck's Law.")
print("B(L, T) = (2 * h * c^2 / L^5) * (1 / (e^(h * c / (L * k * T)) - 1))")
print("\nFirst, let's calculate the theoretical spectral radiance for the given temperature and wavelength.")

# Step 3: Calculate the theoretical spectral radiance based on T and L
exp_numerator = h * c
exp_denominator = L * k * T_measured
exponent = exp_numerator / exp_denominator
B_theoretical = (2 * h * c**2 / L**5) * (1 / (math.exp(exponent) - 1))

# Print the detailed calculation
print(f"\n--- Calculating Theoretical B ---")
print(f"Given L = {L_nm} nm ({L:.2e} m) and T = {T_measured} K:")
print(f"The exponent (h * c / (L * k * T)) is calculated as:")
print(f"({h:.6e} * {c:.6e}) / ({L:.2e} * {k:.6e} * {T_measured}) = {exponent:.4f}")
print(f"Plugging the numbers into Planck's Law:")
print(f"B = (2 * {h:.6e} * {c:.6e}**2 / {L:.2e}**5) * (1 / (e**{exponent:.4f} - 1))")
print(f"Result: The theoretical spectral radiance should be {B_theoretical:.3e} W·sr⁻¹·m⁻²·m⁻¹.")
print(f"This theoretical value is much lower than your measured value of {B_measured:.3e}.\n")

# Step 4: Analyze the discrepancy - can L be the sole error?
print("--- Analyzing the Error ---")
# For a given temperature, there is a maximum possible radiance at a peak wavelength (Wien's Law).
b = 2.898e-3 # Wien's displacement constant in m·K
L_peak = b / T_measured
B_peak = (2 * h * c**2 / L_peak**5) * (1 / (math.exp(h * c / (L_peak * k * T_measured)) - 1))
print(f"At {T_measured} K, the maximum possible radiance is {B_peak:.3e} (at {L_peak*1e9:.0f} nm).")
print(f"Your measured radiance ({B_measured:.3e}) is higher than this maximum.")
print("This proves that the wavelength (L) cannot be the only error. The error must be in T or B.\n")

# Step 5: Calculate the temperature required to produce the measured radiance
print("--- Calculating the Correct Temperature ---")
print("Let's assume your spectral radiance measurement is correct and calculate the required temperature.")
log_argument = (2 * h * c**2 / (L**5 * B_measured)) + 1
T_correct = (h * c / (L * k)) / math.log(log_argument)

# Print the detailed calculation for the corrected temperature
print(f"To get B = {B_measured:.2e} at L = {L:.2e}, the temperature T would need to be:")
print(f"T = (h * c / (L * k)) / ln( (2 * h * c^2 / (L^5 * B)) + 1 )")
print(f"ln_argument = (2 * {h:.6e} * {c:.6e}**2 / ({L:.2e}**5 * {B_measured:.2e})) + 1 = {log_argument:.4f}")
print(f"T = ({h:.6e} * {c:.6e} / ({L:.2e} * {k:.6e})) / ln({log_argument:.4f})")
print(f"Result: The required temperature is {T_correct:.0f} K.")

# Step 6: Use astrophysical context to make a final conclusion
print("\n--- Final Conclusion ---")
print("The star is a DB-class white dwarf. The typical temperature for this class is 12,000 K to 30,000 K.")
print(f"Your estimated temperature of {T_measured} K is outside this range.")
print(f"The corrected temperature of {T_correct:.0f} K is perfectly within the expected range.")
print("\nTherefore, the temperature (T) is the quantity most likely to have a measurement error.")
print(f"The corrected value is {round(T_correct)} K.")

final_answer = f"T{round(T_correct)}"
#The final answer is presented below as requested
#print(f"\nFinal Answer Code: {final_answer}")