# Mapping of group names to their exponents.
group_exponents = {
    "PSL(2,4)": 30,
    "Z4xZ4": 4,
    "D8": 8,
    "S4": 12,
    "A4": 6,
    "D3": 6,
    "Z3xZ3": 3,
    "Z2^3": 2,
}

# Mapping of each visualization (V_i) to its corresponding group name.
# This mapping is derived from analyzing the number of vertices and the graph structure.
graph_to_group = {
    # Column 1
    "V1": "PSL(2,4)",
    "V5": "Z4xZ4",
    "V9": "S4",
    "V13": "Z2^3",
    # Column 2
    "V2": "D8",
    "V6": "Z3xZ3",
    "V10": "Z3xZ3",
    "V14": "PSL(2,4)",
    # Column 3
    "V3": "A4",
    "V7": "D3",
    "V11": "A4",
    "V15": "D8",
    # Column 4
    "V4": "S4",
    "V8": "Z2^3",
    "V12": "D3",
    "V16": "Z4xZ4",
}

# Get the exponents for each graph
exponents = {v: group_exponents[g] for v, g in graph_to_group.items()}

# Define the columns
c1_exps = [exponents["V1"], exponents["V5"], exponents["V9"], exponents["V13"]]
c2_exps = [exponents["V2"], exponents["V6"], exponents["V10"], exponents["V14"]]
c3_exps = [exponents["V3"], exponents["V7"], exponents["V11"], exponents["V15"]]
c4_exps = [exponents["V4"], exponents["V8"], exponents["V12"], exponents["V16"]]

# Calculate the sums
s1 = sum(c1_exps)
s2 = sum(c2_exps)
s3 = sum(c3_exps)
s4 = sum(c4_exps)

# Print the results
print("Calculation of column sums:")
print(f"S1 = {c1_exps[0]} + {c1_exps[1]} + {c1_exps[2]} + {c1_exps[3]} = {s1}")
print(f"S2 = {c2_exps[0]} + {c2_exps[1]} + {c2_exps[2]} + {c2_exps[3]} = {s2}")
print(f"S3 = {c3_exps[0]} + {c3_exps[1]} + {c3_exps[2]} + {c3_exps[3]} = {s3}")
print(f"S4 = {c4_exps[0]} + {c4_exps[1]} + {c4_exps[2]} + {c4_exps[3]} = {s4}")
print("\nThe four column sums as an ordered list are:")
print(f"{{{s1}, {s2}, {s3}, {s4}}}")

# Final answer in the required format
final_answer = f"{{{s1}, {s2}, {s3}, {s4}}}"
# print(f"<<<{final_answer}>>>")