def analyze_patterns(inputs, outputs):
    print("Zero preservation analysis:")
    for i in range(len(inputs)):
        input_list = [int(x) for x in inputs[i].split()]
        output_list = [int(x) for x in outputs[i].split()]
        
        print(f"\nExample {i+1}:")
        print("Input:", input_list)
        print("Output:", output_list)
        
        # Analyze zero preservation
        for j in range(len(input_list)):
            if input_list[j] == 0:
                print(f"Position {j}: Input 0 -> Output {output_list[j]}")
            else:
                print(f"Position {j}: Input {input_list[j]} -> Output {output_list[j]}")

# Example inputs and outputs
inputs = [
    "0 0 5 0 8 8 0 4 0 8",
    "0 6 2 6 0 7 2 0 3 0",
    "9 5 7 3 3 6 0 0 0 3"
]

outputs = [
    "0 0 5 0 5 5 0 5 0 5",
    "0 5 5 5 0 5 5 0 5 0",
    "5 5 5 5 5 5 0 0 0 5"
]

analyze_patterns(inputs, outputs)