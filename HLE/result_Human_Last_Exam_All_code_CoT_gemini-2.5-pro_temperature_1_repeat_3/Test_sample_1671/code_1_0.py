import textwrap

# The goal is to find the NCCN-recommended treatment for HER2+ MBC that has progressed after a regimen containing trastuzumab and a taxane.

# Step 1: Analyze the provided text for the first-line treatment recommendation.
# The text says: "Current National Comprehensive Cancer Network (NCCN) guidelines with category 1 evidence recommend THP as the preferred first-line regimen for HER2+ MBC."
# THP includes a taxane, which matches the condition in the question.

# Step 2: Find the recommended treatment after progression on the first-line therapy.
# The text explicitly states: "Upon disease progression, TDXd is recommended in the second-line setting."

# Step 3: Identify what TDXd stands for.
# The text clarifies: "fam-trastuzumab deruxtecan (TDXd)".

# Step 4: Compare this information with the answer choices.
# A. Tucatinib trastuzumab and capecitabine: The text identifies this as a preferred *third-line* regimen.
# B. Trastuzumab emtansine: The text states this was *replaced* by TDXd as the second-line therapy.
# C. Fam-trastuzumab deruxtecan: This matches the finding from the text as the recommended second-line therapy.
# D. Chemotherapy: This is too general; the guidelines recommend specific agents.

# Step 5: Conclude the correct answer.
conclusion = "Based on the text, Fam-trastuzumab deruxtecan (TDXd) is the NCCN-recommended treatment in the second-line setting after disease progression on a first-line regimen that includes a taxane and trastuzumab."

print(textwrap.fill(conclusion, 80))
print("\nTherefore, the correct choice is C.")
print("<<<C>>>")