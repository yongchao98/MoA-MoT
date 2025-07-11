import math

def analyze_plots():
    """
    Analyzes six plots of quantum evolution to find the physically valid one
    by checking fundamental quantum mechanical constraints.
    """
    print("Analyzing plots based on physical constraints for a single qubit system.")
    print("The three key conditions for a valid evolution are:")
    print("1. -1 <= <σz> <= 1")
    print("2. 0 <= S <= log(2)")
    print("3. <σz>^2 + 4 * |<σ+>|^2 <= 1")
    print("\nChecking a representative point from each plot:\n")

    # Data points extracted from graphs: (plot_name, time, <σz>, |<σ+>|, S)
    # Points are chosen where violations are most likely to be visible.
    test_data = {
        'A': ('A', 2.0, 0.3, 0.9, 0.2),
        'B': ('B', 0.0, 0.5, 0.7, 0.0),
        'C': ('C', 2.0, 1.7, 0.6, -0.8),
        'D': ('D', 5.0, 0.4, 0.4, 0.8),
        'E': ('E', 1.0, 0.72, 0.5, 0.05),
        'F': ('F', 1.5, 0.7, 0.3, 0.15)
    }

    valid_plot = None

    for plot_name in sorted(test_data.keys()):
        name, t, sz, sp_mag, s_val = test_data[plot_name]
        is_valid = True
        
        print(f"--- Checking Plot {name} (at t={t}) ---")

        # Condition 1: Check <σz>
        cond1 = -1 <= sz <= 1
        if not cond1:
            print(f"VIOLATION: <σz> = {sz} is outside the allowed range [-1, 1].")
            is_valid = False

        # Condition 2: Check Entropy S
        cond2_non_negative = s_val >= 0
        max_entropy = math.log(2)
        cond2_upper_bound = s_val <= max_entropy
        if not cond2_non_negative:
            print(f"VIOLATION: Entropy S = {s_val} is negative.")
            is_valid = False
        if not cond2_upper_bound:
            print(f"VIOLATION: Entropy S = {s_val} exceeds the maximum for a qubit (~{max_entropy:.3f}).")
            is_valid = False

        # Condition 3: Check Bloch vector length
        bloch_len_sq = sz**2 + 4 * sp_mag**2
        cond3 = bloch_len_sq <= 1.001  # Tolerance for reading error

        print(f"Checking Bloch vector constraint: <σz>^2 + 4 * |<σ+>|^2 <= 1")
        print(f"  Equation with values: ({sz})^2 + 4 * ({sp_mag})^2")
        print(f"  Calculation: {sz**2:.3f} + 4 * {sp_mag**2:.3f} = {sz**2:.3f} + {4*sp_mag**2:.3f} = {bloch_len_sq:.3f}")
        
        if not cond3:
            print(f"VIOLATION: The result {bloch_len_sq:.3f} is greater than 1.")
            is_valid = False

        if is_valid:
            print("Conclusion: No violations found. This plot appears to be physically valid.")
            valid_plot = name
        else:
            print("Conclusion: Invalid physical evolution.")
        
        print("-"*(21+len(name)) + "\n")

    print("--- Final Result ---")
    if valid_plot:
        print(f"The only plot that represents a physically valid quantum evolution is Plot {valid_plot}.")
        print(f"\n{valid_plot}.")
    else:
        # This case shouldn't be reached as F is valid.
        print("No valid plot was identified.")

if __name__ == '__main__':
    analyze_plots()