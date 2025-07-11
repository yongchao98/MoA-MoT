import numpy as np

def analyze_plots():
    """
    Analyzes the physical validity of the six quantum evolution plots.
    """
    plots = {
        'A': {'sz_range': (-0.2, 0.8), 'sp_abs_range': (0.05, 0.9), 'S_range': (0, 0.45)},
        'B': {'sz_range': (0.5, 0.7), 'sp_abs_range': (0.55, 0.7), 'S_range': (0, 0.25)},
        'C': {'sz_range': (1.0, 1.7), 'sp_abs_range': (0.4, 0.7), 'S_range': (-1.2, 0.4)},
        'D': {'sz_range': (0.35, 0.5), 'sp_abs_range': (0.3, 0.7), 'S_range': (0, 0.8)},
        'E': {'sz_range': (0.5, 0.7), 'sp_abs_range': (0.5, 0.7), 'S_range': (0, 0.2)},
        'F': {'sz_range': (0.5, 0.7), 'sp_abs_range': (0.3, 0.4), 'S_range': (0, 0.25)}
    }

    max_entropy_ln = np.log(2) # approx 0.693

    print("Analyzing the physical constraints for each plot:")
    print("1. Expectation value <σz> must be in [-1, 1].")
    print("2. |<σ+>| must be in [0, 0.5].")
    print(f"3. Entropy S must be in [0, ln(2) ≈ {max_entropy_ln:.3f}].")
    print("-" * 30)

    valid_plot = None
    for name, ranges in plots.items():
        is_valid = True
        reasons = []

        # Check <σz> range
        if ranges['sz_range'][0] < -1 or ranges['sz_range'][1] > 1:
            is_valid = False
            reasons.append(f"<σz> out of range [-1, 1] (max value ~{ranges['sz_range'][1]}).")
        
        # Check |<σ+>| range
        if ranges['sp_abs_range'][0] < 0 or ranges['sp_abs_range'][1] > 0.5:
            is_valid = False
            reasons.append(f"|<σ+>| out of range [0, 0.5] (max value ~{ranges['sp_abs_range'][1]}).")
            
        # Check Entropy S range
        if ranges['S_range'][0] < 0:
            is_valid = False
            reasons.append(f"Entropy S is negative (min value ~{ranges['S_range'][0]}).")
        if ranges['S_range'][1] > max_entropy_ln:
            is_valid = False
            reasons.append(f"Entropy S exceeds max value ln(2) (max value ~{ranges['S_range'][1]}).")

        print(f"Plot {name}:")
        if is_valid:
            print("  - All values are within physical bounds. This plot is physically valid.")
            valid_plot = name
        else:
            print("  - Invalid for the following reasons:")
            for reason in reasons:
                print(f"    - {reason}")
        print()

    if valid_plot:
        print(f"Conclusion: The only physically valid diagram is {valid_plot}.")
    else:
        print("Conclusion: None of the plots appear to be strictly valid based on the estimated ranges.")

analyze_plots()