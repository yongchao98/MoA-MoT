def solve_sequence():
    """
    Solves the sequence puzzle by decomposing it into blocks and finding the pattern.
    """
    # The given sequence
    sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # The proposed next 4 elements from choice F
    next_elements = [2, 3, 2, 1]

    # The full sequence including the proposed answer
    full_sequence = sequence + next_elements

    # Decomposing the full sequence into logical blocks
    blocks = [
        [3, 2, 1],
        [2, 3],
        [3],
        [3, 2],
        [2],
        [2, 3],  # Proposed next block
        [2, 1]   # Proposed final block
    ]
    
    reconstructed_sequence = [num for block in blocks for num in block]

    print("Original sequence: 3 2 1 2 3 3 3 2 2")
    print("The pattern can be understood by breaking the sequence into blocks:")
    
    block_definitions = [
        "B1=(3,2,1): A decreasing sequence from 3.",
        "B2=(2,3): An increasing sequence from 2.",
        "B3=(3): A constant sequence.",
        "B4=(3,2): A decreasing sequence from 3.",
        "B5=(2): A constant sequence.",
        "This suggests the next blocks follow a similar logic.",
        "Based on the pattern, the next block B6 should start at 2 (where B5 ended).",
        "The proposed continuation from choice F is (2,3,2,1). This can be split into two new blocks:",
        "B6=(2,3): An increasing sequence starting at 2.",
        "B7=(2,1): Since B6 ends at boundary 3, B7 starts at 2 and decreases."
    ]
    for definition in block_definitions:
        print(f"- {definition}")
        
    print("\nThis creates a consistent pattern. Therefore, the next 4 elements are:")
    # Print the equation part by part
    print(f"{' '.join(map(str, sequence))} ... -> {' '.join(map(str, next_elements))}")

solve_sequence()
<<<F>>>