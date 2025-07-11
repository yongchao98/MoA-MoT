import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# This script deduces the parameter change for each plot and constructs the final answer string.
# Each print statement explains the reasoning for a single plot's code.

# Plot 1: Baseline
plot1_code = '0'
print(f"Analysis for Plot 1: The plot displays a characteristic asymmetric double-scroll chaotic attractor. This appears to be a standard representation of the system's dynamics and will serve as our baseline or 'initial simulation'.")
print(f"Assigned Code: {plot1_code}")

# Plot 2: Doubled Rayleigh number (R)
plot2_code = 'R'
print(f"\nAnalysis for Plot 2: Compared to Plot 1, the attractor here is significantly denser and the time-series oscillations are of a higher frequency. This suggests more vigorous chaotic dynamics, which is a classic effect of increasing the Rayleigh number (R), the primary driving parameter. This corresponds to a doubled Rayleigh number.")
print(f"Assigned Code: {plot2_code}")

# Plot 3: Changed initial condition (z)
plot3_code = 'z'
print(f"\nAnalysis for Plot 3: The overall structure of the attractor is virtually identical to that of Plot 1. For a chaotic system, starting from a different initial condition (like Z0) leads to a different path along the same attractor. This explains the strong similarity. This corresponds to a change in Z0.")
print(f"Assigned Code: {plot3_code}")

# Plot 4: Halved mu (m)
plot4_code = 'm'
print(f"\nAnalysis for Plot 4: A striking feature of this plot is the high degree of symmetry between the blue and orange attractors, which contrasts with the asymmetry in Plot 1. The parameter μ governs the symmetry of the coupling terms. Bringing μ closer to 1 increases symmetry. This effect matches the result of halving μ.")
print(f"Assigned Code: {plot4_code}")

# Plot 5: Halved Rayleigh number (r)
plot5_code = 'r'
print(f"\nAnalysis for Plot 5: The dynamics have fundamentally changed from chaotic to stable. The trajectories are shown spiraling into a fixed point, confirmed by the damped oscillations in the time series. This bifurcation is caused by reducing the Rayleigh number below the chaotic threshold. This corresponds to halving R.")
print(f"Assigned Code: {plot5_code}")

# Plot 6: Halved Prandtl number (p)
plot6_code = 'p'
print(f"\nAnalysis for Plot 6: The attractor is simpler and sparser than the baseline, and the time-series oscillations are noticeably slower. The Prandtl number (P) controls the time scale of the z and Z variables. Halving P would make these variables 'slower', simplifying the overall dynamics as seen here.")
print(f"Assigned Code: {plot6_code}")


# Combine the codes into the final string
final_string = plot1_code + plot2_code + plot3_code + plot4_code + plot5_code + plot6_code
print(f"\nCombining the codes for plots 1 through 6 in order gives the final string: {final_string}")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the captured output
print(output)
print("<<<0Rzmrp>>>")