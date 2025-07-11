import collections

# The mapping was determined by ranking the patterns in each group by their
# visual density (number of black cells) and pairing them based on that rank.
#
# Density Ranking for Patterns A-H (least to most dense):
# D, G, A, B, E, F, C, H
#
# Density Ranking for Patterns 1-8 (least to most dense):
# 3, 4, 1, 6, 8, 7, 2, 5
#
# This leads to the following one-to-one mapping:
# D -> 3
# G -> 4
# A -> 1
# B -> 6
# E -> 8
# F -> 7
# C -> 2
# H -> 5

# Create a dictionary to hold the mapping
mapping = {
    'A': 1,
    'B': 6,
    'C': 2,
    'D': 3,
    'E': 8,
    'F': 7,
    'G': 4,
    'H': 5
}

# The problem asks for the result in a specific order: A, B, C, D, E, F, G, H
ordered_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
result_list = [mapping[label] for label in ordered_labels]

# Format the output as requested: {N_A,N_B,N_C,N_D,N_E,N_F,N_G,N_H}
result_string = "{" + ",".join(map(str, result_list)) + "}"

print(f"The established mapping is A->{mapping['A']}, B->{mapping['B']}, C->{mapping['C']}, D->{mapping['D']}, E->{mapping['E']}, F->{mapping['F']}, G->{mapping['G']}, H->{mapping['H']}.")
print("The final answer in the required format is:")
print(result_string)

# Final answer in the specified format for the platform
final_answer_for_submission = result_string
# print(f"<<<{final_answer_for_submission}>>>")