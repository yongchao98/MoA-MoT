import sys
# Redirect print to a string to control the final output format
from io import StringIO
old_stdout = sys.stdout
sys.stdout = mystdout = StringIO()


def solve_protein_puzzle():
    """
    Analyzes a series of protein interaction experiments to determine the most correct conclusion.
    """
    # Theoretical masses (kDa) from amino acid sequences
    theo_mass = {'A': 25, 'B': 150, 'C': 60, 'D': 100, 'Kinase': 40}

    # --- Step 1: Analyze Experiment 1 (Individual Proteins) ---
    print("--- Analysis of Experiment 1 (Individual Proteins) ---")
    obs_mass_exp1 = {'A': 50, 'B': 300, 'C': 60, 'D': 100}

    # Calculate and print the oligomeric state for each protein
    a_state = obs_mass_exp1['A'] / theo_mass['A']
    b_state = obs_mass_exp1['B'] / theo_mass['B']
    c_state = obs_mass_exp1['C'] / theo_mass['C']
    d_state = obs_mass_exp1['D'] / theo_mass['D']

    print(f"Protein A: Observed mass {obs_mass_exp1['A']} kDa / Theoretical mass {theo_mass['A']} kDa = {a_state:.0f}. State: Dimer.")
    print(f"Protein B: Observed mass {obs_mass_exp1['B']} kDa / Theoretical mass {theo_mass['B']} kDa = {b_state:.0f}. State: Dimer.")
    print(f"Protein C: Observed mass {obs_mass_exp1['C']} kDa / Theoretical mass {theo_mass['C']} kDa = {c_state:.0f}. State: Monomer.")
    print(f"Protein D: Observed mass {obs_mass_exp1['D']} kDa / Theoretical mass {theo_mass['D']} kDa = {d_state:.0f}. State: Monomer.")
    print("Conclusion from Exp 1: Protein A and B are stable dimers, while C and D are monomers.\n")

    # Define masses for use in subsequent steps
    a_dimer_mass = obs_mass_exp1['A']
    b_dimer_mass = obs_mass_exp1['B']
    c_mono_mass = obs_mass_exp1['C']
    d_mono_mass = obs_mass_exp1['D']

    # --- Step 2: Analyze Experiment 2 (All Proteins Mixed) ---
    print("--- Analysis of Experiment 2 (Proteins Mixed without Kinase) ---")
    peaks_exp2 = [300, 210]
    print(f"Observed peaks: {peaks_exp2[0]} kDa and {peaks_exp2[1]} kDa.")
    print(f"The {peaks_exp2[0]} kDa peak corresponds to unbound Protein B dimer ({b_dimer_mass} kDa).")
    
    complex_mass_exp2 = peaks_exp2[1]
    acd_complex_calc = a_dimer_mass + c_mono_mass + d_mono_mass
    print(f"Let's test the composition of the {complex_mass_exp2} kDa peak.")
    print(f"Equation: Protein A (dimer) + Protein C (monomer) + Protein D (monomer) = {a_dimer_mass} kDa + {c_mono_mass} kDa + {d_mono_mass} kDa = {acd_complex_calc} kDa.")
    print("The calculated mass matches the observed peak.")
    print("Conclusion from Exp 2: In a direct competition, non-phosphorylated Protein A has a higher affinity for C+D than Protein B does.\n")

    # --- Step 3: Analyze Experiment 3 (Mixed with Kinase) ---
    print("--- Analysis of Experiment 3 (Mixed with Kinase) ---")
    peaks_exp3 = [25, 40, 460]
    print(f"Observed peaks: {peaks_exp3[0]} kDa, {peaks_exp3[1]} kDa, and {peaks_exp3[2]} kDa.")
    print(f"The {peaks_exp3[1]} kDa peak is the Kinase ({theo_mass['Kinase']} kDa).")
    print(f"The {peaks_exp3[0]} kDa peak is Protein A monomer ({theo_mass['A']} kDa), indicating phosphorylation caused its dimer to dissociate.")

    complex_mass_exp3 = peaks_exp3[2]
    bcd_complex_calc = b_dimer_mass + c_mono_mass + d_mono_mass
    print(f"Let's test the composition of the {complex_mass_exp3} kDa peak.")
    print(f"Equation: Protein B (dimer) + Protein C (monomer) + Protein D (monomer) = {b_dimer_mass} kDa + {c_mono_mass} kDa + {d_mono_mass} kDa = {bcd_complex_calc} kDa.")
    print("The calculated mass matches the observed peak.")
    print("Conclusion from Exp 3: Phosphorylation of Protein A decreases its affinity for C+D, allowing Protein B to bind instead.\n")

    # --- Step 4: Analyze Experiment 4 (Phosphatase Treatment) ---
    print("--- Analysis of Experiment 4 (Dephosphorylation) ---")
    peaks_exp4 = [50, 460]
    print(f"Observed peaks: {peaks_exp4[0]} kDa and {peaks_exp4[1]} kDa.")
    print(f"The {peaks_exp4[0]} kDa peak is Protein A ({a_dimer_mass} kDa) re-forming its dimer after dephosphorylation.")
    print(f"The {peaks_exp4[1]} kDa peak shows the Protein B+C+D complex ({bcd_complex_calc} kDa) remains intact.")
    print("Conclusion from Exp 4: The B+C+D complex is very stable, as dephosphorylated Protein A cannot displace Protein B.\n")
    
    # --- Step 5: Final Evaluation of Answer Choices ---
    print("--- Final Conclusion ---")
    print("Based on the combined evidence:")
    print("1. From Exp 2, non-phosphorylated Protein A has a higher binding affinity for C+D than Protein B does.")
    print("2. From all experiments, Protein B is consistently a dimer (300 kDa).")
    print("3. Protein A's oligomeric state changes (dimer or monomer), while C and D are monomers.")
    print("These facts directly support answer choice G: 'Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A, protein B always exists as a dimer in opposition to proteins A, C, and D.'")
    
# Run the analysis
solve_protein_puzzle()

# Restore stdout and print the captured output along with the final answer format
sys.stdout = old_stdout
print(mystdout.getvalue())
print("<<<G>>>")