import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Create a new string IO object
sys.stdout = captured_output = io.StringIO()

# The logic for determining the six-character string is as follows:
# We designate one plot as the reference ('0') and identify the single parameter
# change that results in each of the other five plots.

# Plot 1: Assumed to be the reference ('0') simulation, showing a typical chaotic attractor.
code_1 = '0'

# Plot 2: Shows a much denser and larger attractor with higher frequency oscillations.
# This indicates increased chaoticity due to a stronger driving force.
# This corresponds to doubling the Rayleigh number (R).
code_2 = 'R'

# Plot 3: The attractor has the same size and shape as Plot 1, but with a different trajectory.
# For a chaotic system, this is the classic signature of a different initial condition (Z0).
code_3 = 'Z'

# Plot 4: This attractor is significantly more symmetric and synchronized than Plot 1.
# The parameter μ controls the asymmetry. Halving μ ('m') would make the system
# more symmetric (assuming base μ > 1), which explains the visual change.
code_4 = 'm'

# Plot 5: The trajectory spirals into a stable fixed point, meaning the chaos has been suppressed.
# This happens when the driving force is reduced below the chaotic threshold.
# This corresponds to halving the Rayleigh number (r).
code_5 = 'r'

# Plot 6: The system has transitioned from chaos to a simpler, periodic orbit.
# In Lorenz-like systems, increasing the Prandtl number (P) often simplifies the dynamics,
# leading to such periodic limit cycles.
code_6 = 'P'

# Assemble the final six-character string.
final_string = code_1 + code_2 + code_3 + code_4 + code_5 + code_6

# Restore original stdout
sys.stdout = original_stdout
# Get the content of captured_output
output = captured_output.getvalue()

# Print the final combined string according to the identified parameters.
# The instruction "output each number in the final equation" is interpreted as
# outputting each character of the final answer string with its justification.
print(f"The identified code for each plot is:")
print(f"Plot 1: Character '{code_1}'. Reason: This plot serves as the reference 'initial simulation' with its characteristic chaotic attractor.")
print(f"Plot 2: Character '{code_2}'. Reason: The Rayleigh number (R) is doubled, increasing the system's driving force and leading to a more vigorous, dense chaotic motion.")
print(f"Plot 3: Character '{code_3}'. Reason: The initial condition (Z0) is changed. On a chaotic attractor, this results in a different trajectory but preserves the attractor's overall shape.")
print(f"Plot 4: Character '{code_4}'. Reason: The temperature ratio (μ) is halved. This change brings the system closer to symmetry, resulting in the observed highly synchronized and symmetric behavior.")
print(f"Plot 5: Character '{code_5}'. Reason: The Rayleigh number (r) is halved, reducing the driving force below the threshold for chaos and causing the system to decay to a stable fixed point.")
print(f"Plot 6: Character '{code_6}'. Reason: The Prandtl number (P) is doubled. This slows down the velocity components relative to thermal effects, simplifying the dynamics into a periodic orbit.")
print(f"\nAssembling the characters in order gives the final string: {final_string}")
print(f"<<<{final_string}>>>")