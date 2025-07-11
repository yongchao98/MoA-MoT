# The analysis for each CFG is done step-by-step above.
# This code block simply prints the final formatted result based on that analysis.

# Profile for X1: [S, s, 33]
# S: Scheme
# s: separated
# 33: dimension

# Profile for X2: [DM, s, irr, 3]
# DM: Deligne-Mumford stack
# s: separated
# irr: irreducible
# 3: dimension

# Profile for X3: [S, s, 7]
# S: Scheme
# s: separated
# 7: dimension

x1_profile = "[S, s, 33]"
x2_profile = "[DM, s, irr, 3]"
x3_profile = "[S, s, 7]"

print(f"{x1_profile} {x2_profile} {x3_profile}")