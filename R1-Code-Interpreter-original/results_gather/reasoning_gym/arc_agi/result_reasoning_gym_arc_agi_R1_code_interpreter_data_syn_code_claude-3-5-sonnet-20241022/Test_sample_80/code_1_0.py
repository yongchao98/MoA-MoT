def find_nine_pattern(examples):
    patterns = []
    for example in examples:
        output = [list(map(int, row.split())) for row in example.split('\n')]
        nine_pos = []
        for i in range(len(output)):
            for j in range(len(output[i])):
                if output[i][j] == 9:
                    nine_pos.append((i,j))
        patterns.append(nine_pos)
    return patterns

# All example outputs
example_outputs = [
"""3 3 3
3 3 3
3 3 9""",

"""3 3 3
3 3 3
3 3 9""",

"""9 3 3
3 9 3
3 3 9""",

"""3 3 3
3 9 3
3 3 9""",

"""3 3 3
3 9 3
3 3 9"""
]

patterns = find_nine_pattern(example_outputs)
print("Nine positions in each example:")
for i, pattern in enumerate(patterns, 1):
    print(f"Example {i}: {pattern}")