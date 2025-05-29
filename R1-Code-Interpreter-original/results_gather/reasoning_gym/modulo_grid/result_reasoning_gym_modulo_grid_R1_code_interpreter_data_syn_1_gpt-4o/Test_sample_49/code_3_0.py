# Calculate the differences between consecutive column indices
column_differences = [cols[i] - cols[i - 1] for i in range(1, len(cols))]

# Print the differences
print("Column Differences:", column_differences)