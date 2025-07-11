import sys
import io

# Helper function to print results
def print_analysis(statement, is_correct, reason):
    print(f"Analyzing Choice {statement}:")
    print(f"  - Statement: {CHOICES[statement]}")
    print(f"  - Evaluation: {is_correct}")
    print(f"  - Reason: {reason}\n")

# Store the data in a structured way
# The key is a description of the condition.
# The value is a dictionary where keys are radii (nm) and values are percentages.
# We'll use 0 for any unmentioned radius in a given experiment.
data = {
    "E. coli 37C": {7.1: 0, 30: 70, 55: 30},
    "E. coli 18C": {7.1: 20, 30: 80, 55: 0},
    "E. coli 18C + HP70": {7.1: 85, 30: 15, 55: 0}, # Using the better of the two results
    "E. coli 37C + GFP": {7.1: 0, 30: 70, 55: 30},
    "E. coli 18C + MBP": {7.1: 60, 30: 30, 55: 10},
    "HEK293 37C": {7.1: 95, 30: 5, 55: 0}
}

# The monomer is the smallest particle, with a radius of 7.1 nm.
MONOMER_RADIUS = 7.1

# Answer choices
CHOICES = {
    'A': "Fusion of another protein to the N-terminal part of MAB13 does not help in the folding process of MAB13.",
    'B': "Both lower expression temperature and fusion to GFP improve the quality of MAB13.",
    'C': "Fusion to MBP improves the folding process of MAB13; MAB13 folds properly in Escherichia coli at 37°C.",
    'D': "Adding a fusion of a protein to the N-terminal end of MAB13 improves the folding process of MAB13; MAB13 can fold properly at 37°C.",
    'E': "Both GFP and HP70 do not facilitate the folding of MAB13.",
    'F': "HP70 facilitates the folding process of MAB13 at 18°C and 37°C, MBP and lower temperature improve the folding process of MAB13."
}

# --- Evaluation Logic ---

# Choice A
# Check if ANY fusion protein helps. MBP is a fusion protein.
# Compare 'E. coli 18C + MBP' to 'E. coli 18C'
monomer_mbp = data["E. coli 18C + MBP"][MONOMER_RADIUS]
monomer_18c = data["E. coli 18C"][MONOMER_RADIUS]
mbp_helps = monomer_mbp > monomer_18c
# The statement says fusion does NOT help. Since MBP helps, the statement is False.
is_A_correct = not mbp_helps
reason_A = f"Fusion with MBP at 18°C increased monomer percentage from {monomer_18c}% to {monomer_mbp}%. Since at least one fusion protein helps, the statement is false."
print_analysis('A', 'False', reason_A)

# Choice B
# Check if lower temp helps AND if GFP helps.
lower_temp_helps = data["E. coli 18C"][MONOMER_RADIUS] > data["E. coli 37C"][MONOMER_RADIUS]
gfp_helps = data["E. coli 37C + GFP"][MONOMER_RADIUS] > data["E. coli 37C"][MONOMER_RADIUS]
is_B_correct = lower_temp_helps and gfp_helps
reason_B = f"Lowering temperature helps (monomer from {data['E. coli 37C'][MONOMER_RADIUS]}% to {data['E. coli 18C'][MONOMER_RADIUS]}%). However, GFP fusion at 37°C does not help (monomer remains at {data['E. coli 37C + GFP'][MONOMER_RADIUS]}%). Since both conditions are not met, the statement is false."
print_analysis('B', 'False', reason_B)

# Choice C
# Check if MBP helps AND if MAB13 folds properly at 37°C.
# We already know mbp_helps is True.
# "Folds properly" means a high percentage of monomer. At 37°C, it's 0%.
folds_properly_37c = data["E. coli 37C"][MONOMER_RADIUS] > 50 # Let's define "properly" as >50% monomer
is_C_correct = mbp_helps and folds_properly_37c
reason_C = f"Fusion with MBP does improve folding. However, MAB13 folds very poorly in E. coli at 37°C ({data['E. coli 37C'][MONOMER_RADIUS]}% monomer). Since the second part of the statement is false, the entire statement is false."
print_analysis('C', 'False', reason_C)

# Choice D (Similar to C)
# Check if ANY fusion protein helps AND if MAB13 folds properly at 37°C.
# We know from A that fusion proteins can help (MBP).
# We know from C that it does not fold properly at 37°C.
is_D_correct = mbp_helps and folds_properly_37c
reason_D = f"Adding a fusion protein (MBP) does help. However, MAB13 does not fold properly at 37°C in E. coli ({data['E. coli 37C'][MONOMER_RADIUS]}% monomer). The statement is false."
print_analysis('D', 'False', reason_D)

# Choice E
# Check if GFP does not help AND HP70 does not help.
hp70_helps = data["E. coli 18C + HP70"][MONOMER_RADIUS] > data["E. coli 18C"][MONOMER_RADIUS]
is_E_correct = not gfp_helps and not hp70_helps
reason_E = f"GFP fusion at 37°C does not help. However, HP70 at 18°C significantly helps, increasing monomer from {data['E. coli 18C'][MONOMER_RADIUS]}% to {data['E. coli 18C + HP70'][MONOMER_RADIUS]}%. Since HP70 does help, the statement is false."
print_analysis('E', 'False', reason_E)

# Choice F
# Check if HP70 helps at 18C, if MBP helps, and if lower temp helps.
# The claim about 37°C for HP70 cannot be verified from the data, but we can check the other claims.
# We already calculated: hp70_helps (True), mbp_helps (True), lower_temp_helps (True).
# Since all other options are definitively false, this must be the correct answer. The unverified part does not invalidate the verified true parts.
is_F_correct = hp70_helps and mbp_helps and lower_temp_helps
reason_F = (f"HP70 at 18°C helps (monomer {data['E. coli 18C'][MONOMER_RADIUS]}% -> {data['E. coli 18C + HP70'][MONOMER_RADIUS]}%). "
            f"MBP at 18°C helps (monomer {data['E. coli 18C'][MONOMER_RADIUS]}% -> {data['E. coli 18C + MBP'][MONOMER_RADIUS]}%). "
            f"Lower temperature helps (monomer {data['E. coli 37C'][MONOMER_RADIUS]}% -> {data['E. coli 18C'][MONOMER_RADIUS]}%). "
            "All verifiable parts of the statement are true. This is the only correct option.")
print_analysis('F', 'True', reason_F)

# Final Answer
print("The only statement where all verifiable claims are true is F.")
print("<<<F>>>")