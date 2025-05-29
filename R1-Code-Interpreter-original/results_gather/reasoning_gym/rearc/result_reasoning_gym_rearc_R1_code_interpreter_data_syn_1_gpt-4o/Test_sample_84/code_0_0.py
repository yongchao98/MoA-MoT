def extract_first_elements(input_grid):
    # Convert input grid to a list of lists
    input_rows = [list(map(int, row.split())) for row in input_grid.strip().split('\n')]
    
    # Use a set to track seen rows
    seen_rows = set()
    output = []
    
    for row in input_rows:
        row_tuple = tuple(row)
        if row_tuple not in seen_rows:
            seen_rows.add(row_tuple)
            output.append(row[0])
    
    return output

# Test input grid
test_input = """
1
1
5
"""

# Get the output grid
output_grid = extract_first_elements(test_input)

# Print the output grid
print("\n".join(map(str, output_grid)))