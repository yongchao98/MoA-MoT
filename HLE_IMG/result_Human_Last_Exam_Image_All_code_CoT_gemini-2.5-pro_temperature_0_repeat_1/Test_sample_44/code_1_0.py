import sys
import io

# Redirect stdout to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_nmr_puzzle():
    """
    This function analyzes 1H NMR data to identify the corresponding compound from a list of choices.
    """
    # Step 1: Analyze the provided 1H NMR data.
    nmr_data = {
        8.19: 1, 7.79: 1, 7.47: 1, 7.38: 1, 6.98: 1, 6.63: 1, 6.61: 1,
        4.19: 4, 3.63: 4, 3.21: 2, 2.83: 2, 1.98: 2
    }
    
    integrations = list(nmr_data.values())
    
    # Calculate total observed protons
    total_observed_protons = sum(integrations)
    
    # Calculate observed aromatic protons (typically > 6.5 ppm)
    aromatic_protons_observed = sum(integration for shift, integration in nmr_data.items() if shift > 6.5)
    
    print("Step 1: Analysis of the 1H NMR Data")
    print("=====================================")
    print("The integrations for the observed signals are: " + ", ".join(map(str, integrations)))
    integration_sum_str = " + ".join(map(str, integrations))
    print(f"The total number of observed protons is the sum of these integrations:")
    print(f"{integration_sum_str} = {total_observed_protons} protons.")
    print(f"The number of observed aromatic protons (with chemical shifts > 6.5 ppm) is {aromatic_protons_observed}.")
    print("\n")

    # Step 2: Analyze the structures of the candidate compounds.
    print("Step 2: Analysis of the Candidate Compounds")
    print("==========================================")
    # Protons in structural fragments
    quinoline_arom_H = 3
    quinoline_aliph_H = 6  # 3 x CH2
    piperazine_H = 8       # 4 x CH2
    thiourea_NH = 1
    pyridine_arom_H = 4
    phenyl_arom_H = 5

    # Calculate theoretical proton counts for each compound
    protons_A = {"total": quinoline_arom_H + quinoline_aliph_H + thiourea_NH + piperazine_H + pyridine_arom_H, "aromatic": quinoline_arom_H + pyridine_arom_H}
    protons_C = {"total": quinoline_arom_H + quinoline_aliph_H + thiourea_NH + piperazine_H + phenyl_arom_H, "aromatic": quinoline_arom_H + phenyl_arom_H}
    protons_B_E = {"total": 2 * protons_A["total"], "aromatic": 2 * protons_A["aromatic"]} # B and E are isomers
    protons_D = {"total": 2 * protons_C["total"], "aromatic": 2 * protons_C["aromatic"]}
    
    compounds = {"A": protons_A, "B": protons_B_E, "C": protons_C, "D": protons_D, "E": protons_B_E}

    for name, data in sorted(compounds.items()):
        print(f"Compound {name}: Expected Total Protons = {data['total']}, Expected Aromatic Protons = {data['aromatic']}")
    print("\n")

    # Step 3: Compare NMR data with compound structures and deduce the answer.
    print("Step 3: Comparison and Deduction")
    print("================================")
    print("Let's compare the observed data with the expected values for each compound.")
    
    print("\n- Compounds B, D, and E are metal complexes with two ligands.")
    print(f"  They have expected total proton counts of {protons_B_E['total']} or {protons_D['total']}, which is much higher than the observed {total_observed_protons} protons.")
    print("  Therefore, we can eliminate B, D, and E.")

    print("\n- Now let's compare the free ligands, A and C.")
    print(f"  Compound A: Has {protons_A['aromatic']} aromatic protons. This MATCHES the {aromatic_protons_observed} aromatic protons observed in the NMR spectrum.")
    print(f"  Compound C: Has {protons_C['aromatic']} aromatic protons. This DOES NOT MATCH the {aromatic_protons_observed} aromatic protons observed.")
    
    print("\n- Let's confirm with the total proton count for Compound A.")
    print(f"  Compound A has a theoretical total of {protons_A['total']} protons.")
    print(f"  The NMR spectrum shows {total_observed_protons} protons.")
    print(f"  The difference is {protons_A['total']} - {total_observed_protons} = 1 proton.")
    print("  This difference is explained by the exchangeable thiourea N-H proton, which is often not observed or is too broad to be integrated in a 1H NMR spectrum.")
    print("\n")

    # Step 4: Final Conclusion
    print("Step 4: Conclusion")
    print("==================")
    print("The NMR data, especially the count of 7 aromatic protons and a total of 21 observed protons, is fully consistent with the structure of Compound A.")
    
    # The answer choices are A. B, B. D, C. E, D. C, E. A
    print("The correct compound is A. Looking at the answer choices, this corresponds to choice E.")

solve_nmr_puzzle()

# Get the captured output and print it to the actual stdout
output = captured_output.getvalue()
sys.stdout = old_stdout
print(output)
print("<<<E>>>")