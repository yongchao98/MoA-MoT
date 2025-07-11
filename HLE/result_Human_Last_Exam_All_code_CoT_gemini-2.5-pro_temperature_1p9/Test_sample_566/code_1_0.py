import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to capture the output
sys.stdout = captured_output = io.StringIO()

# --- Start of the Explanation Block ---
# In Coherent Anti-Stokes Raman Scattering (CARS), the key frequency relationship
# for generating the signal is:
# ω_as = ω_p + ω_pr - ω_s
# where ω_as is the anti-Stokes frequency, ω_p is the pump, ω_pr is the probe, and ω_s is the Stokes.

# The signal is resonantly enhanced when the beat frequency between the pump and Stokes
# matches a molecular vibration frequency (Ω):
# Ω = ω_p - ω_s

# In broadband CARS, one of the beams (typically the Stokes, but can be the pump)
# has a wide frequency bandwidth. This allows for the simultaneous excitation of a
# whole range of vibrational frequencies (Ω).

# Consequently, the generated anti-Stokes signal (ω_as) is also broadband.
# When this signal is analyzed with a spectrometer, its spectrum (intensity vs. frequency)
# reveals peaks corresponding to the different molecular vibrations that were excited.

# Therefore, the broadband anti-Stokes beam contains distinct, resolvable spectral features
# that provide information about the chemical composition of the sample.

print("CARS Signal Equation: anti_stokes = pump + probe - stokes")
print("CARS Resonance Condition: pump - stokes = vibrational_frequency")

final_answer_choice = 'C'
final_answer_text = "You can generate a anti Stokes beam that contains distinguishable information"

print(f"\nBased on the principles of broadband CARS, the correct statement is:")
print(f"Choice: {final_answer_choice}")
print(f"Text: {final_answer_text}")

# --- End of the Explanation Block ---

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the final captured output
print(output)