import collections

# Step 1: Define group assignments and their exponents
# The key is the visualization number (1-16)
# The value is a tuple: (Group Name, Exponent)
group_data = {
    1: ("PSL(2,4)", 30),
    2: ("Z4 x Z4", 4),
    3: ("A4", 6),
    4: ("S4", 12),
    5: ("Z2^3", 2),
    6: ("D8", 8),
    7: ("D3", 6),
    8: ("Z3 x Z3", 3),
    9: ("PSL(2,4)", 30),
    10: ("Z3 x Z3", 3),
    11: ("A4", 6),
    12: ("D3", 6),
    13: ("Z4 x Z4", 4),
    14: ("S4", 12),
    15: ("Z2^3", 2),
    16: ("D8", 8)
}

# Step 2: Define the grid layout
# Columns are 1-indexed for clarity
grid_columns = {
    1: [1, 5, 9, 13],
    2: [2, 6, 10, 14],
    3: [3, 7, 11, 15],
    4: [4, 8, 12, 16]
}

# Step 3: Calculate and print the sum for each column
column_sums = []
print("Calculating the sum of exponents for each column:\n")

for col_num in sorted(grid_columns.keys()):
    v_indices = grid_columns[col_num]
    
    # Get the exponents for the current column
    exponents = [group_data[v_idx][1] for v_idx in v_indices]
    
    # Calculate the sum
    total_sum = sum(exponents)
    column_sums.append(total_sum)
    
    # Format the equation string
    equation_str = " + ".join(map(str, exponents))
    
    # Print the result for the current column
    print(f"S_{col_num} = exp(V_{v_indices[0]}) + exp(V_{v_indices[1]}) + exp(V_{v_indices[2]}) + exp(V_{v_indices[3]})")
    print(f"S_{col_num} = {equation_str} = {total_sum}\n")

# Step 4: Print the final ordered list of sums
print("The four column sums as an ordered list are:")
print(column_sums)
