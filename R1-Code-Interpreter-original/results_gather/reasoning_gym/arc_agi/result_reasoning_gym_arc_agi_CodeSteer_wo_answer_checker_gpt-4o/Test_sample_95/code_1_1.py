def extract_rules(input_examples, output_examples):
    rules = []
    for input_grid, output_grid in zip(input_examples, output_examples):
        grid_rules = []
        for input_row, output_row in zip(input_grid, output_grid):
            row_rules = []
            for i, (input_val, output_val) in enumerate(zip(input_row, output_row)):
                if input_val == output_val:
                    row_rules.append(('constant', output_val))
                elif output_val == min(input_row):
                    row_rules.append(('min', i))
                elif output_val == max(input_row):
                    row_rules.append(('max', i))
                else:
                    row_rules.append(('unknown', None))
            grid_rules.append(row_rules)
        rules.append(grid_rules)
    return rules

def apply_rules(input_grid, rules):
    output_grid = []
    for row, row_rules in zip(input_grid, rules[0]):  # Using the first set of rules for simplicity
        new_row = []
        for i, (rule_type, value) in enumerate(row_rules):
            if rule_type == 'constant':
                new_row.append(value)
            elif rule_type == 'min':
                new_row.append(min(row))
            elif rule_type == 'max':
                new_row.append(max(row))
            else:
                new_row.append(row[i])  # Default to input value if rule is unknown
        output_grid.append(new_row)
    return output_grid

# Example input-output pairs
input_examples = [
    [[2, 7, 7], [2, 7, 7], [9, 7, 7], [7, 7, 7], [7, 0, 9], [7, 0, 7], [7, 9, 7], [7, 7, 7], [9, 4, 7], [7, 4, 7], [7, 4, 7]],
    [[7, 5, 7], [7, 5, 7], [5, 5, 7], [9, 7, 7], [7, 7, 7], [7, 7, 9], [7, 7, 6], [7, 7, 9], [7, 7, 7], [9, 4, 7], [7, 4, 7]],
    [[7, 4, 7], [7, 4, 7], [7, 9, 7], [7, 7, 7], [9, 6, 7], [7, 9, 7], [7, 7, 7], [7, 4, 9], [7, 4, 7]],
    [[7, 7, 7], [7, 7, 4], [7, 4, 4], [7, 9, 7], [7, 7, 7], [9, 4, 7], [7, 9, 7], [7, 7, 7], [7, 7, 9], [7, 6, 6], [7, 6, 7]]
]

output_examples = [
    [[7, 2, 7], [7, 2, 7], [7, 2, 7], [0, 0, 7], [0, 7, 7], [0, 7, 7], [4, 4, 7], [7, 4, 7], [7, 4, 7]],
    [[7, 5, 7], [7, 5, 7], [5, 5, 7], [5, 7, 7], [6, 7, 7], [6, 7, 7], [6, 7, 7], [4, 4, 7], [7, 4, 7]],
    [[7, 4, 7], [7, 4, 7], [7, 4, 7], [7, 6, 6], [7, 7, 6], [7, 4, 4], [7, 4, 7]],
    [[7, 7, 4], [7, 4, 4], [7, 4, 7], [7, 4, 4], [7, 7, 4], [7, 7, 6], [7, 6, 6], [7, 6, 7]]
]

# Extract rules from examples
rules = extract_rules(input_examples, output_examples)

# Test input grid
test_input_grid = [
    [7, 0, 0],
    [7, 7, 9],
    [7, 7, 7],
    [7, 5, 9],
    [7, 9, 7],
    [7, 7, 7],
    [9, 6, 6],
    [7, 7, 9],
    [7, 7, 7],
    [7, 4, 9],
    [7, 4, 7]
]

# Apply rules to the test input grid
predicted_output_grid = apply_rules(test_input_grid, rules)

# Print the predicted output grid
for row in predicted_output_grid:
    print(' '.join(map(str, row)))