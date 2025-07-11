# The analysis leads to the following mapping from condition to simulation index:
# Condition 1 (minimum t) is met by Simulation 2.
# Condition 2 (minimum |s|) is met by Simulation 1.
# Condition 3 (unique sign(s)) is met by Simulation 4.
# Condition 4 (maximum s) is met by Simulation 4.

# The required output is the sequence of these simulation indices.
c1 = 2
c2 = 1
c3 = 4
c4 = 4

# Print the result as a single number as requested by the format.
# "Answer with the simulation indices ordered by the condition met (e.g., ... 4321)"
print(f"{c1}{c2}{c3}{c4}")
