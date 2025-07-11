import sys
import io

# This problem requires analyzing the provided images based on knowledge of dynamical systems,
# specifically Lorenz-like equations. The code's purpose is to output the final answer string
# derived from this analysis.

def solve_lorenz_plots():
    """
    Deduces the six-character string corresponding to the parameter changes in plots 1-6.

    The reasoning is as follows:
    1. Plot 1: Appears to be the baseline chaotic attractor. Code: '0'.
    2. Plot 2: The attractor is much denser and the time-series oscillations are of higher frequency.
               This indicates a stronger driving force, characteristic of an increased Rayleigh number. Code: 'R'.
    3. Plot 3: The attractor is almost identical in shape and character to Plot 1. This is the expected result
               of changing only the initial conditions on a chaotic attractor. Code: 'z' (or 'Z').
    4. Plot 4: The system's behavior has fundamentally changed from chaotic to periodic and synchronized.
               The time series shows perfectly regular oscillations. This bifurcation is often caused by
               strong coupling, which is controlled by the Biot number. Code: 'B'.
    5. Plot 5: The attractor exhibits very tight spirals around the two lobes. This happens when the z-dynamics
               are slow relative to the x-y dynamics, which corresponds to a large Prandtl number. Code: 'P'.
    6. Plot 6: The attractor is "sparser" and less complex than the baseline in Plot 1. This suggests
               weaker chaotic behavior, consistent with a reduced driving force from a halved Rayleigh number. Code: 'r'.

    Combining these deductions gives the sequence for plots 1 through 6.
    """
    
    # The final string is constructed by concatenating the codes for each plot in order.
    plot_1_code = '0'  # Baseline
    plot_2_code = 'R'  # Rayleigh doubled
    plot_3_code = 'z'  # Z0 halved (initial condition)
    plot_4_code = 'B'  # Biot doubled
    plot_5_code = 'P'  # Prandtl doubled
    plot_6_code = 'r'  # Rayleigh halved

    final_string = plot_1_code + plot_2_code + plot_3_code + plot_4_code + plot_5_code + plot_6_code
    
    print(final_string)

solve_lorenz_plots()