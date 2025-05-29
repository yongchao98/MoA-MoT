# Calculate the differences between consecutive column positions
differences = [columns[i+1] - columns[i] for i in range(len(columns) - 1)]

# Print the differences to analyze the pattern
print("Differences:", differences)