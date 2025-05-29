def combine_blocks():
    # Initial counts of blocks
    blocks = {
        '[A]': 2,
        '[B]': 7,
        '[C]': 2,
        '{A}': 0,
        '{B}': 0,
        '{C}': 0,
        '(A)': 0,
        '(B)': 0,
        '(C)': 0
    }

    # Define transformation rules
    transformations = [
        ({'[A]': 1, '[B]': 1, '[C]': 1}, {'{A}': 1}),
        ({'[A]': 1, '[B]': 1}, {'{C}': 1}),
        ({'[B]': 1, '[C]': 1}, {'{B}': 1}),
        ({'[C]': 2}, {'{C}': 1}),
        ({'{A}': 1, '{C}': 1}, {'(A)': 1, '(B)': 1}),
        ({'{B}': 2}, {'(C)': 1})
    ]

    # Function to apply transformations
    def apply_transformations():
        nonlocal blocks
        changes = False

        for inputs, outputs in transformations:
            # Check if the transformation can be applied
            if all(blocks[block] >= count for block, count in inputs.items()):
                # Apply the transformation
                for block, count in inputs.items():
                    blocks[block] -= count
                for block, count in outputs.items():
                    blocks[block] += count
                changes = True

        return changes

    # Apply transformations until no more changes can be made
    while apply_transformations():
        pass

    # Output the final counts
    result = ', '.join(f"{count} {block}" for block, count in blocks.items())
    print(result)

combine_blocks()