import math

def solve_edx_problem():
    """
    Determines the lightest element detectable in a specific EDX setup by considering
    detector window transmission and spectral overlaps from the sample matrix.
    """

    # --- Setup Parameters ---
    beam_energy_keV = 30.0
    
    # Be Window properties
    window_thickness_um = 100.0
    window_thickness_cm = window_thickness_um / 10000.0  # convert µm to cm
    be_density_g_cm3 = 1.85

    # Sample matrix properties
    matrix_element = "Tungsten (W)"
    w_interference_lines = {
        "W Mα": 1.774  # keV
    }

    # Detector property
    detector_resolution_keV = 0.130 # A typical detector energy resolution

    # Candidate elements from choices
    # Data: Kα energy (keV), K-edge energy (keV), Mass Attenuation Coeff of Be at Kα energy (cm²/g)
    # Mass attenuation data from Lawrence Berkeley National Lab (LBL) CXRO database.
    elements = {
        "Na": {"Z": 11, "Ka_keV": 1.041, "K_edge_keV": 1.072, "mu_Be_cm2_g": 1618},
        "Mg": {"Z": 12, "Ka_keV": 1.254, "K_edge_keV": 1.303, "mu_Be_cm2_g": 966},
        "Si": {"Z": 14, "Ka_keV": 1.740, "K_edge_keV": 1.839, "mu_Be_cm2_g": 423},
        "Ca": {"Z": 20, "Ka_keV": 3.691, "K_edge_keV": 4.038, "mu_Be_cm2_g": 65.6}
    }
    
    print("Analyzing EDX detectability for light elements on a Tungsten (W) sample.")
    print(f"Setup: {beam_energy_keV} keV electron beam, SSD with {window_thickness_um} µm Be window.\n")
    
    detectable_elements = []

    # Sort elements by atomic number to evaluate from lightest to heaviest
    sorted_elements = sorted(elements.items(), key=lambda item: item[1]['Z'])

    for name, properties in sorted_elements:
        print(f"--- Analyzing: {name} (Z={properties['Z']}) ---")
        
        # --- 1. Excitation Check ---
        can_be_excited = properties['K_edge_keV'] < beam_energy_keV
        print(f"Step 1: Excitation Check")
        if not can_be_excited:
            print(f"  Result: CANNOT be excited ({properties['K_edge_keV']} keV K-edge > {beam_energy_keV} keV beam). Skipping {name}.\n")
            continue
        else:
            print(f"  Result: Can be excited ({properties['K_edge_keV']} keV K-edge < {beam_energy_keV} keV beam).")

        # --- 2. Transmission Check ---
        mu = properties['mu_Be_cm2_g']
        exponent = -mu * be_density_g_cm3 * window_thickness_cm
        transmission = math.exp(exponent)
        
        print(f"Step 2: Transmission Check through {window_thickness_um} µm Be window")
        print(f"  Transmission Formula: T = exp(-μ * ρ * x)")
        print(f"  {name} Kα Energy = {properties['Ka_keV']:.3f} keV")
        print(f"  μ (Be Mass Attenuation Coeff at {properties['Ka_keV']:.3f} keV) = {mu} cm²/g")
        print(f"  ρ (Be Density) = {be_density_g_cm3} g/cm³")
        print(f"  x (Window Thickness) = {window_thickness_cm} cm")
        print(f"  Calculation: T = exp(-{mu} * {be_density_g_cm3} * {window_thickness_cm}) = exp({exponent:.3f})")
        print(f"  Result: Transmission = {transmission:.3e} ({transmission:.4%})")

        transmission_threshold = 0.0001 # A practical threshold of 0.01%
        if not transmission > transmission_threshold:
            print(f"  Conclusion: Transmission is negligible. {name} is NOT detectable.\n")
            continue
        else:
            print(f"  Conclusion: Transmission is sufficient for potential detection.")

        # --- 3. Interference Check ---
        print(f"Step 3: Spectral Interference Check with {matrix_element} matrix")
        is_distinguishable = True
        for line_name, line_energy in w_interference_lines.items():
            energy_diff_keV = abs(properties['Ka_keV'] - line_energy)
            print(f"  Checking {name} Kα ({properties['Ka_keV']:.3f} keV) against {line_name} ({line_energy:.3f} keV).")
            print(f"  Equation: |E_elem - E_matrix| = |{properties['Ka_keV']:.3f} - {line_energy:.3f}| = {energy_diff_keV:.3f} keV")
            print(f"  Detector Resolution ≈ {detector_resolution_keV:.3f} keV")
            if energy_diff_keV < detector_resolution_keV:
                is_distinguishable = False
                break
        
        if not is_distinguishable:
            print(f"  Conclusion: Peaks overlap. {name} is NOT distinguishable from the W matrix.\n")
            continue
        else:
             print(f"  Conclusion: Peaks are well-separated. {name} is distinguishable.\n")

        # --- Conclusion for this element ---
        detectable_elements.append(name)
        print(f"Overall Conclusion for {name}: DETECTABLE.\n")


    # --- Final Conclusion ---
    print("--- FINAL RESULT ---")
    if not detectable_elements:
        print("Based on this analysis, none of the candidate light elements are detectable.")
    else:
        lightest_element = min(detectable_elements, key=lambda x: elements[x]['Z'])
        print(f"The only element from the list that passes all checks is: {lightest_element}")
        
solve_edx_problem()
<<<D>>>