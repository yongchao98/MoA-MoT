import pandas as pd

def analyze_relaxation_data():
    """
    Analyzes the effect of methylation on ring relaxation time based on data from the plots.
    """
    # Estimated data for ring 1 from the plots
    data = {
        'Temperature (K)': [325, 340, 350, 365, 400],
        'N1_tau_ring1 (ns)': [200, 50, 25, 10, 3],  # Nonmethylated
        'M1_tau_ring1 (ns)': [150, 30, 12, 8, 2]   # Methylated
    }
    
    df = pd.DataFrame(data)
    
    print("--- Analysis of Relaxation Times (Part 1) ---")
    print("Comparing the relaxation time <τ> of ring 1 for Nonmethylated (N1) vs. Methylated (M1) systems.\n")
    
    for index, row in df.iterrows():
        temp = row['Temperature (K)']
        n1_tau = row['N1_tau_ring1 (ns)']
        m1_tau = row['M1_tau_ring1 (ns)']
        comparison = "lower" if m1_tau < n1_tau else "higher"
        print(f"At {temp} K: τ(N1) ≈ {n1_tau} ns, τ(M1) ≈ {m1_tau} ns. The methylated ring has a {comparison} relaxation time.")

    print("\nConclusion for Part 1: The data shows that the addition of a methyl group decreases the relaxation time, meaning the ring rotates faster.")
    
    print("\n--- Analysis of Transition Temperature (Part 2) ---")
    print("The addition of a bulky lateral methyl group disrupts the ability of molecules to pack efficiently.")
    print("This steric hindrance weakens the intermolecular forces that stabilize the ordered nematic phase.")
    print("Conclusion for Part 2: As a result, the nematic phase becomes less stable, and the nematic-isotropic transition temperature is expected to decrease.")
    
    print("\n--- Final Evaluation ---")
    print("Option E correctly combines these two points:")
    print("1. Decreased relaxation time due to steric effects (supported by the data).")
    print("2. Lower nematic-isotropic transition temperature due to disrupted molecular packing.")

analyze_relaxation_data()