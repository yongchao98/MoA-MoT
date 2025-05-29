def analyze_changes(input_str, output_str):
    input_arr = [int(x) for x in input_str.split()]
    output_arr = [int(x) for x in output_str.split()]
    
    # Print positions where changes occur
    changes = []
    for i in range(len(input_arr)):
        if input_arr[i] != output_arr[i]:
            changes.append(i)
    
    # Analyze sequences of 2's and their transformations
    sequences = []
    current_seq = []
    for i in range(len(input_arr)):
        if input_arr[i] == 2:
            current_seq.append(i)
        elif current_seq:
            sequences.append(current_seq)
            current_seq = []
    if current_seq:
        sequences.append(current_seq)
    
    print("Sequences of 2's positions:", sequences)
    print("Changed positions:", changes)

# Analyze all examples
examples = [
    ("0 0 2 2 2 2 0 2 2 2 0 0 0 0 0 0 0 2 2 2 0 0 0 0 0 0 0 0",
     "0 0 8 8 2 2 0 8 2 2 0 0 0 0 0 0 0 8 2 2 0 0 0 0 0 0 0 0"),
    ("0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 2 2 0 0 0 0 0 0 0 0 0",
     "0 8 8 8 8 8 8 8 2 2 2 2 2 2 2 0 0 8 2 0 0 0 0 0 0 0 0 0"),
    ("2 2 2 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 2 2 0 0 0 0 0 0 0",
     "8 2 2 0 8 8 8 8 8 8 8 2 2 2 2 2 2 2 0 8 2 0 0 0 0 0 0 0")
]

for i, (inp, out) in enumerate(examples):
    print(f"\nExample {i+1}:")
    analyze_changes(inp, out)