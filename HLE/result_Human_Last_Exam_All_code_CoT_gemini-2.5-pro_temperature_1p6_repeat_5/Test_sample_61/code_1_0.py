# The user wants a comma-separated list of the relevant sections and subsections
# of the UK's Employment Rights Act 1996.

# 1. The section that grants the power to make the regulation: 43FA(3)
# 2. The general section governing regulations: 236
# 3. The subsection confirming the use of a statutory instrument: 236(1)
# 4. The subsection listing powers for the affirmative procedure (which must be checked and ruled out): 236(3)
# 5. The subsection that applies the negative procedure by default: 236(4)

relevant_sections = ["43FA(3)", "236(1)", "236(3)", "236(4)"]

# Print the final comma-separated list as requested.
print(','.join(relevant_sections))
