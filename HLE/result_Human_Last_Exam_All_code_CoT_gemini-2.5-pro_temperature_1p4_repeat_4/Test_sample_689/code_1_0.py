import collections

def analyze_mab13_folding():
    """
    Analyzes DLS data for MAB13 folding and evaluates multiple-choice answers.
    """
    # Store the experimental data
    # Structure: {condition: [(radius, intensity), ...]}
    data = {
        "E. coli 37°C": [(30, 70), (55, 30)],
        "E. coli 18°C": [(7.1, 20), (30, 80)],
        "E. coli 18°C + HP70": [(7.1, 70), (30, 30)], # Note: second HP70 result (85%) confirms the trend
        "HEK293 37°C (Control)": [(7.1, 95), (30, 5)],
        "E. coli 37°C + GFP": [(30, 70), (55, 30)],
        "E. coli 18°C + MBP": [(7.1, 60), (30, 30), (55, 10)]
    }

    MONOMER_RH = 7.1

    # Step 1: Explain the benchmark for good folding
    print("--- Analysis of Protein Folding ---")
    print(f"Based on the control experiment in HEK293 cells, a properly folded MAB13 monomer has a hydrodynamic radius of {MONOMER_RH} nm.")
    print("A higher intensity percentage at this radius indicates better folding and less aggregation.\n")

    # Helper function to get monomer percentage
    def get_monomer_percent(results):
        for radius, intensity in results:
            if radius == MONOMER_RH:
                return intensity
        return 0

    # Step 2: Create a summary table of monomer percentages
    print("--- Summary of Correctly Folded Monomer ---")
    monomer_percentages = collections.OrderedDict()
    for condition, results in data.items():
        monomer_percentages[condition] = get_monomer_percent(results)
    
    for condition, percentage in monomer_percentages.items():
        print(f"{condition:<25}: {percentage}%")
    print("\n--- Evaluation of Answer Choices ---")

    # Step 3: Evaluate each choice
    # A
    ecoli_18c = monomer_percentages["E. coli 18°C"]
    mbp_18c = monomer_percentages["E. coli 18°C + MBP"]
    print("A. Fusion of another protein to the N-terminal part of MAB13 does not help...")
    print(f"   - FALSE. MBP fusion at 18°C increased monomer from {ecoli_18c}% to {mbp_18c}%.\n")

    # B
    ecoli_37c = monomer_percentages["E. coli 37°C"]
    gfp_37c = monomer_percentages["E. coli 37°C + GFP"]
    print("B. Both lower expression temperature and fusion to GFP improve the quality...")
    print(f"   - FALSE. Lowering temperature helped ({ecoli_37c}% -> {ecoli_18c}%), but GFP fusion at 37°C did not ({gfp_37c}% monomer).\n")
    
    # C
    print("C. Fusion to MBP improves the folding...; MAB13 folds properly in E. coli at 37°C.")
    print(f"   - FALSE. MAB13 folds poorly in E. coli at 37°C ({ecoli_37c}% monomer).\n")

    # D
    hek_37c = monomer_percentages["HEK293 37°C (Control)"]
    print("D. Adding a fusion... improves folding...; MAB13 can fold properly at 37°C.")
    print(f"   - FALSE. This is misleading. It folds at 37°C in HEK293 cells ({hek_37c}%), not in E. coli ({ecoli_37c}%).\n")

    # E
    hp70_18c = monomer_percentages["E. coli 18°C + HP70"]
    print("E. Both GFP and HP70 do not facilitate the folding...")
    print(f"   - FALSE. HP70 clearly facilitates folding, increasing monomer from {ecoli_18c}% to {hp70_18c}%.\n")

    # F
    print("F. HP70 facilitates... at 18°C and 37°C, MBP and lower temperature improve...")
    print("   - Let's check the claims based on the data:")
    print(f"     1. HP70 facilitates at 18°C: TRUE ({ecoli_18c}% -> {hp70_18c}%)")
    print("     2. Data for HP70 at 37°C: Not Provided.")
    print(f"     3. Lower temperature improves folding: TRUE ({ecoli_37c}% -> {ecoli_18c}%)")
    print(f"     4. MBP improves folding: TRUE ({ecoli_18c}% -> {mbp_18c}%)")
    print("   - CONCLUSION: This is the only statement where the verifiable claims are all true. All other statements are demonstrably false.")

analyze_mab13_folding()