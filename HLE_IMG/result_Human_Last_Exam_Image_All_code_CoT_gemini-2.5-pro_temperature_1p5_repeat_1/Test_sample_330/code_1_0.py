import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# This script deduces the modified parameter for each of the six plots.
# The logic is based on a visual analysis of the dynamical systems presented.

# Plot 1: Assumed to be the baseline reference simulation. Code '0'.
plot1_code = '0'

# Plot 2: Shows much more chaotic and dense attractors compared to Plot 1.
# This increased chaoticity is characteristic of a higher Rayleigh number (R). Code 'R'.
plot2_code = 'R'

# Plot 3: The attractor shape is very similar to Plot 1, but the density
# is shifted heavily to the orange system. This suggests a different trajectory on the
# same underlying attractor, which is caused by a change in initial conditions (Z0). Code 'Z'.
plot3_code = 'Z'

# Plot 4: Exhibits strong synchronization and symmetry between the two systems.
# This is caused by increasing the coupling strength, controlled by the Biot number (B). Code 'B'.
plot4_code = 'B'

# Plot 5: The system's oscillations decay, and it settles into a stable fixed point.
# This transition from chaos to a stable state is caused by decreasing the
# Rayleigh number (R) below the chaotic threshold. Code 'r'.
plot5_code = 'r'

# Plot 6: The attractors are less "loopy" and more "skeletal" than in Plot 1.
# This structural change, with fewer spirals, is characteristic of a lower
# Prandtl number (P), which speeds up the z-dynamics. Code 'p'.
plot6_code = 'p'

# The final answer is a concatenation of the codes for plots 1 through 6.
final_answer_string = plot1_code + plot2_code + plot3_code + plot4_code + plot5_code + plot6_code

print(final_answer_string)

# Restore stdout
sys.stdout = old_stdout
# Print the captured output to the actual console
output = captured_output.getvalue()
print(output)
print(f'<<<{output.strip()}>>>')