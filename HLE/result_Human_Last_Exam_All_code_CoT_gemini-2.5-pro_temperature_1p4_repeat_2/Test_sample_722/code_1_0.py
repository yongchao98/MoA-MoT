def analyze_sec_mals_data():
    """
    Analyzes the provided SEC-MALS data to determine protein interactions
    and evaluates the given answer choices.
    """
    # Theoretical masses from amino acid sequence
    theoretical_mass = {
        'A': 25,
        'B': 150,
        'C': 60,
        'D': 100
    }
    kinase_mass = 40

    print("### Step-by-Step Analysis ###\n")

    # --- Experiment 1: Analysis of Individual Proteins ---
    print("--- Analysis of Experiment 1: Individual Proteins ---")
    exp1_mass = {
        'A': 50,
        'B': 300,
        'C': 60,
        'D': 100
    }
    
    # Determine oligomeric states
    A_dimer = exp1_mass['A']
    B_dimer = exp1_mass['B']
    C_mono = exp1_mass['C']
    D_mono = exp1_mass['D']

    print(f"Protein A: Theoretical mass={theoretical_mass['A']} kDa, Experimental mass={exp1_mass['A']} kDa. This means Protein A forms a dimer ({exp1_mass['A']} / {theoretical_mass['A']} = 2).")
    print(f"Protein B: Theoretical mass={theoretical_mass['B']} kDa, Experimental mass={exp1_mass['B']} kDa. This means Protein B forms a dimer ({exp1_mass['B']} / {theoretical_mass['B']} = 2).")
    print(f"Protein C: Theoretical mass={theoretical_mass['C']} kDa, Experimental mass={exp1_mass['C']} kDa. This means Protein C is a monomer ({exp1_mass['C']} / {theoretical_mass['C']} = 1).")
    print(f"Protein D: Theoretical mass={theoretical_mass['D']} kDa, Experimental mass={exp1_mass['D']} kDa. This means Protein D is a monomer ({exp1_mass['D']} / {theoretical_mass['D']} = 1).")
    print("\nConclusion from Exp 1: A and B are dimers; C and D are monomers.\n")

    # --- Experiment 2: Mixture of A, B, C, D ---
    print("--- Analysis of Experiment 2: Mixture of A, B, C, D ---")
    exp2_peaks = [300, 210]
    print(f"Observed peaks at {exp2_peaks[0]} kDa and {exp2_peaks[1]} kDa.")
    
    # Identify complexes
    peak1_id = "Protein B Dimer"
    complex_2_mass = A_dimer + C_mono + D_mono
    print(f"The 300 kDa peak corresponds to the unbound Protein B dimer.")
    print(f"Let's test combinations for the 210 kDa peak:")
    print(f"  Protein A (Dimer) + Protein C (Monomer) + Protein D (Monomer) = {A_dimer} kDa + {C_mono} kDa + {D_mono} kDa = {complex_2_mass} kDa.")
    print("\nConclusion from Exp 2: A complex of (A+C+D) forms. Protein B does not participate, suggesting that non-phosphorylated Protein A has a higher affinity for the C+D pair than Protein B does under these initial conditions.\n")

    # --- Experiment 3: Mixture with Kinase ---
    print("--- Analysis of Experiment 3: Mixture with Kinase ---")
    exp3_peaks = [25, 40, 460]
    print(f"Observed peaks at {exp3_peaks[0]} kDa, {exp3_peaks[1]} kDa, and {exp3_peaks[2]} kDa.")
    
    # Identify components and complexes
    A_mono_phos = theoretical_mass['A']
    complex_3_mass = B_dimer + C_mono + D_mono
    print(f"The 40 kDa peak is the unbound Kinase.")
    print(f"The 25 kDa peak corresponds to the theoretical mass of monomeric Protein A. This means the kinase phosphorylated Protein A, causing its dimer to break apart.")
    print(f"Let's test combinations for the 460 kDa peak:")
    print(f"  Protein B (Dimer) + Protein C (Monomer) + Protein D (Monomer) = {B_dimer} kDa + {C_mono} kDa + {D_mono} kDa = {complex_3_mass} kDa.")
    print("\nConclusion from Exp 3: Phosphorylation of Protein A causes it to become a monomer and dissociate from Proteins C and D. This allows Protein B to bind to C and D, forming a new complex.\n")
    
    # --- Experiment 4: Phosphatase Treatment ---
    print("--- Analysis of Experiment 4: Phosphatase Treatment ---")
    exp4_peaks = [50, 460]
    print(f"Observed peaks at {exp4_peaks[0]} kDa and {exp4_peaks[1]} kDa.")
    
    # Identify components
    print(f"The sample from Exp 3 was treated to remove the phosphate from Protein A.")
    print(f"The 50 kDa peak corresponds to the mass of the Protein A dimer. This confirms that dephosphorylation causes Protein A to re-form its dimer.")
    print(f"The 460 kDa peak is the (B+C+D) complex, which remains stable and intact.")
    print("\nConclusion from Exp 4: Dephosphorylated Protein A dimerizes again. However, it cannot displace Protein B from the very stable (B+C+D) complex once it has formed.\n")
    
    # --- Final Evaluation ---
    print("--- Final Evaluation of Answer Choices ---")
    print("Based on the analysis:")
    print("A, C, H are incorrect because the evidence points to phosphorylation of Protein A, not B.")
    print("B and D are incorrect because phosphorylation of Protein A *decreases* its affinity for C and D, it does not increase it.")
    print("E is incorrect because it claims Protein A *always* exists as a homodimer, but it exists as a monomer when phosphorylated.")
    print("F is incorrect because it claims Protein A exists as a monomer (implying always), which is only true when phosphorylated.")
    print("G is incorrect because it claims Protein B *never* has a higher affinity than nonphosphorylated Protein A, which is contradicted by the stability of the B+C+D complex in Experiment 4.")
    print("Since choices A through I all contain at least one statement that is demonstrably false based on the experimental data, none of them can be the correct answer.")
    print("\nTherefore, the only correct option is J.")

analyze_sec_mals_data()