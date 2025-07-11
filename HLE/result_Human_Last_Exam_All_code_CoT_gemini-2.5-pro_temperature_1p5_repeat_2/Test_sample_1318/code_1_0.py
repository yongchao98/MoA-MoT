# Data from Experiment 1
# cfu/ul (plaque-forming units per microliter)
cfu_no_rp_wt = 100000
cfu_no_rp_deltaXY = 100000
cfu_with_rp_wt = 80000
cfu_with_rp_deltaXY = 40000

# --- Calculations ---

# 1. Calculate the resistance provided by the RP system.
# A factor > 1 means the RP system reduces the phage's effectiveness.
resistance_factor_on_wt = cfu_no_rp_wt / cfu_with_rp_wt
resistance_factor_on_deltaXY = cfu_no_rp_deltaXY / cfu_with_rp_deltaXY

# 2. Calculate the anti-defense benefit provided by the XY operon.
# This measures how much better the wt phage is than the deltaXY phage in the presence of the RP system.
# A factor > 1 means the XY operon helps the phage overcome the defense system.
anti_defense_factor_xy = cfu_with_rp_wt / cfu_with_rp_deltaXY

# --- Output Results and Justification ---

print("--- Analysis of Experiment 1 Data ---")
print("This analysis supports the first part of statement H: 'System RP increases the resistance of the bacteria against phageDE3.'")
print(f"Resistance effect of RP on wt phage: {cfu_no_rp_wt} / {cfu_with_rp_wt} = {resistance_factor_on_wt:.2f}-fold reduction in plaques.")
print(f"Resistance effect of RP on deltaXY phage: {cfu_no_rp_deltaXY} / {cfu_with_rp_deltaXY} = {resistance_factor_on_deltaXY:.2f}-fold reduction in plaques.")
print("\nSince the plaque counts are lower in the presence of RP, the RP system increases bacterial resistance.")

print(f"\nAnti-defense effect of XY operon: {cfu_with_rp_wt} / {cfu_with_rp_deltaXY} = {anti_defense_factor_xy:.2f}-fold increase in plaques.")
print("This shows the XY operon helps the phage counteract the RP system's defense.")

print("\n--- Analysis of Experiment 2 Data ---")
print("This analysis supports the second part of statement H: '...because the enzymes XY1 or XY2 can synthesize their products only in the presence of the system RP.'")
print("The molecule of 500 Da (the product) was only detected in Sample 1, which had both:")
print("1. The phage with the XY operon (PhageDE3-wt)")
print("2. The bacteria with the RP system")
print("The absence of either component resulted in no product being detected. This shows the synthesis is dependent on the presence of the RP system.")

print("\n--- Conclusion ---")
print("Statement H correctly combines the key findings from both experiments into a single coherent explanation of the biological interaction.")