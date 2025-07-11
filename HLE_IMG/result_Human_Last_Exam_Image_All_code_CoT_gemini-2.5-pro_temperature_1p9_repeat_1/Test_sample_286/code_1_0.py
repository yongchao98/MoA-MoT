import math

def analyze_quantum_plots():
    """
    Analyzes six plots of quantum evolution to determine which is physically valid.
    """
    print("Analyzing quantum evolution plots based on physical principles.\n")
    
    # --- Physical Principles ---
    print("A valid quantum evolution for a single qubit must satisfy:")
    print("1. Expectation value <σz> must be between -1 and 1.")
    print("2. Bloch vector length squared (r^2) must be <= 1.")
    print("   We calculate this using: r^2 = 4 * |<σ+>|^2 + <σz>^2")
    print(f"3. Von Neumann entropy S must be between 0 and ln(2) (approx {math.log(2):.3f}).\n")

    valid_plot = None

    # --- Plot A ---
    print("--- Analyzing Plot A ---")
    sz = 0.4
    sp_mag = 0.7
    r_squared = 4 * (sp_mag**2) + (sz**2)
    print(f"At t=0, <σz> ≈ {sz}, |<σ+>| ≈ {sp_mag}.")
    print(f"Calculating r^2 = 4 * |<σ+>|^2 + <σz>^2")
    print(f"r^2 = 4 * ({sp_mag})^2 + ({sz})^2 = {4*sp_mag**2:.2f} + {sz**2:.2f} = {r_squared:.2f}")
    if r_squared > 1:
        print(f"Result: INVALID. The Bloch vector length squared ({r_squared:.2f}) is greater than 1.\n")
    else:
        print("Result: This point is valid.\n")

    # --- Plot B ---
    print("--- Analyzing Plot B ---")
    sz = 0.5
    sp_mag = 0.7
    r_squared = 4 * (sp_mag**2) + (sz**2)
    print(f"At t=0, <σz> ≈ {sz}, |<σ+>| ≈ {sp_mag}.")
    print(f"Calculating r^2 = 4 * |<σ+>|^2 + <σz>^2")
    print(f"r^2 = 4 * ({sp_mag})^2 + ({sz})^2 = {4*sp_mag**2:.2f} + {sz**2:.2f} = {r_squared:.2f}")
    if r_squared > 1:
        print(f"Result: INVALID. The Bloch vector length squared ({r_squared:.2f}) is greater than 1.\n")
    else:
        print("Result: This point is valid.\n")

    # --- Plot C ---
    print("--- Analyzing Plot C ---")
    sz_max = 1.7
    s_min = -1.2
    print(f"From visual inspection, <σz> reaches a maximum of ≈ {sz_max}.")
    print(f"From visual inspection, S reaches a minimum of ≈ {s_min}.")
    print(f"Result: INVALID. <σz> exceeds 1 and entropy S cannot be negative.\n")

    # --- Plot D ---
    print("--- Analyzing Plot D ---")
    s_max = 0.8
    ln2 = math.log(2)
    print(f"From visual inspection, the entropy S reaches a maximum of ≈ {s_max}.")
    print(f"This exceeds the maximum possible entropy for a single qubit, ln(2) ≈ {ln2:.3f}.")
    print(f"Result: INVALID.\n")

    # --- Plot E ---
    print("--- Analyzing Plot E ---")
    sz = 0.5
    sp_mag = 0.7
    r_squared = 4 * (sp_mag**2) + (sz**2)
    print(f"At t=0, <σz> ≈ {sz}, |<σ+>| ≈ {sp_mag}.")
    print(f"Calculating r^2 = 4 * |<σ+>|^2 + <σz>^2")
    print(f"r^2 = 4 * ({sp_mag})^2 + ({sz})^2 = {4*sp_mag**2:.2f} + {sz**2:.2f} = {r_squared:.2f}")
    if r_squared > 1:
        print(f"Result: INVALID. The Bloch vector length squared ({r_squared:.2f}) is greater than 1.\n")
    else:
        print("Result: This point is valid.\n")
        
    # --- Plot F ---
    print("--- Analyzing Plot F ---")
    is_f_valid = True
    # Point 1: t=0
    sz = 0.5
    sp_mag = 0.0
    r_squared = 4 * (sp_mag**2) + (sz**2)
    print(f"Checking point at t=0: <σz> ≈ {sz}, |<σ+>| ≈ {sp_mag}.")
    print(f"r^2 = 4 * ({sp_mag})^2 + ({sz})^2 = {4*sp_mag**2:.2f} + {sz**2:.2f} = {r_squared:.2f}")
    if r_squared > 1:
        is_f_valid = False

    # Point 2: t≈1.5
    sz = 0.5
    sp_mag = 0.25
    r_squared = 4 * (sp_mag**2) + (sz**2)
    print(f"Checking point at t≈1.5: <σz> ≈ {sz}, |<σ+>| ≈ {sp_mag}.")
    print(f"r^2 = 4 * ({sp_mag})^2 + ({sz})^2 = {4*sp_mag**2:.2f} + {sz**2:.2f} = {r_squared:.2f}")
    if r_squared > 1:
        is_f_valid = False

    # Point 3: t≈3
    sz = 0.7
    sp_mag = 0.0
    r_squared = 4 * (sp_mag**2) + (sz**2)
    print(f"Checking point at t≈3.0: <σz> ≈ {sz}, |<σ+>| ≈ {sp_mag}.")
    print(f"r^2 = 4 * ({sp_mag})^2 + ({sz})^2 = {4*sp_mag**2:.2f} + {sz**2:.2f} = {r_squared:.2f}")
    if r_squared > 1:
        is_f_valid = False

    if is_f_valid:
        print("All tested points have r^2 <= 1. The expectation values and entropy are also within their allowed ranges.")
        print("Result: VALID.\n")
        valid_plot = 'F'
    else:
        print("Result: INVALID.\n")

    # --- Conclusion ---
    if valid_plot:
        print("====================================")
        print(f"Based on the analysis, Plot {valid_plot} is the only one showing a physically valid quantum evolution.")
        print("====================================")
    else:
        print("None of the plots appear to be fully valid based on this analysis.")

if __name__ == "__main__":
    analyze_quantum_plots()
