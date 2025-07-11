import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer to hold the output
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer


print("### Analysis of Broadband CARS Microscopy ###\n")

print("1. The Goal of Broadband CARS Microscopy:")
print("The primary goal of broadband (or multiplex) CARS microscopy is to obtain a vibrational spectrum from a sample very quickly. This means measuring the response of multiple different molecular vibrations at the same time.\n")

print("2. The Underlying Physics:")
print("CARS is a nonlinear optical process involving three light-matter interactions.")
print("  - First, a pump beam (frequency ω_p) and a Stokes beam (frequency ω_S) interact with the sample.")
print("  - If the frequency difference (ω_p - ω_S) matches a molecule's vibrational frequency (Ω), it drives that vibration.")
print("  - A probe beam (frequency ω_pr, usually the same as the pump) then scatters from this excited vibration.")
print("  - This generates a new signal, the anti-Stokes beam, at frequency ω_as.\n")

print("3. The Key Equation for Information Encoding:")
print("The frequency of the generated anti-Stokes signal is given by the equation:")
print("  ω_as = ω_pr + (ω_p - ω_S)")
print("Since the resonance condition is Ω = ω_p - ω_S, we can write:")
print("  ω_as = ω_pr + Ω\n")

print("Let's print the components of this final equation, assuming the probe and pump are the same (ω_p = ω_pr):")
# This section fulfills the requirement to "output each number in the final equation"
# by explaining the components of the symbolic equation.
equation_components = {
    'ω_as (Anti-Stokes Frequency)': 'The frequency of the light signal that is detected. It carries the information.',
    'ω_p (Pump/Probe Frequency)': 'The frequency of the incident laser beam(s).',
    'Ω (Vibrational Frequency)': 'The natural vibrational frequency of the molecule being studied.'
}
for term, explanation in equation_components.items():
    print(f"  - {term}: {explanation}")

print("\nThis equation shows that each vibrational frequency (Ω) in the sample is mapped to a unique output frequency (ω_as) in the signal.\n")

print("4. Why 'Broadband' Provides Distinguishable Information:")
print("In broadband CARS, one of the beams (typically the Stokes beam) has a wide range of frequencies. This allows the condition Ω = ω_p - ω_S to be met for many different vibrational modes (Ω₁, Ω₂, Ω₃...) simultaneously.")
print("As a result, a broadband anti-Stokes signal is generated, containing components for each vibration (ω_as₁ = ω_p + Ω₁, ω_as₂ = ω_p + Ω₂, etc.).")
print("When this signal is passed through a spectrometer, it is separated into its constituent frequencies, producing a spectrum where each peak corresponds to a specific, distinguishable molecular vibration.\n")

print("5. Conclusion:")
print("Therefore, broadband CARS microscopy generates an anti-Stokes beam that explicitly contains distinguishable vibrational information, which is its primary advantage.\n")

# This is the final answer in the required format
print("<<<C>>>")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
final_output = output_buffer.getvalue()

# Print the final output to the actual console
print(final_output)