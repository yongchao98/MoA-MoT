def predict_output(inp):
    # Find sequences
    def get_sequences(arr):
        result = []
        current = arr[0]
        count = 1
        start = 0
        for i in range(1, len(arr)):
            if arr[i] != current:
                result.append((start, current, count))
                current = arr[i]
                count = 1
                start = i
            else:
                count += 1
        result.append((start, current, count))
        return result
    
    sequences = get_sequences(inp)
    print("Input sequences:", sequences)
    
    # Check if all numbers are same
    if len(sequences) == 1:
        return inp
    
    # Check if it's a single zero in sequence of same numbers
    if len(sequences) == 3 and sequences[1][1] == 0 and sequences[1][2] == 1:
        if sequences[0][1] == sequences[2][1]:  # Same number before and after zero
            return [sequences[0][1]] * 4 + [0] + [sequences[0][1]] * 19
    
    # Check if it's zeros surrounded by same numbers
    if len(sequences) == 3 and sequences[1][1] == 0:
        if sequences[0][1] == sequences[2][1]:  # Same number before and after zeros
            return [0] * 5 + [sequences[0][1]] * 17 + [0] * 2
    
    # If no pattern matches, return input unchanged
    return inp

test_input = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,6,0,0,0,0,0]
print("\nPredicted output:", predict_output(test_input))