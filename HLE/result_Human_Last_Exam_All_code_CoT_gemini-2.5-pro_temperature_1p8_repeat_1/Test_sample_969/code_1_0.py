# The plan is to programmatically generate the sequence based on a discovered pattern.
# The pattern consists of two parts: a sequence of values and a sequence of counts.
# 1. The values 'v' bounce between 3 and 1 (e.g., 3, 2, 1, 2, 3, 2, ...).
# 2. The count 'c' for how many times a value is repeated is 1 for the first four values,
#    and 3 for all subsequent values.
# The given sequence is an incomplete part of this generated sequence.

def solve_sequence():
    """
    This function generates the sequence based on the derived rules
    and prints the next 4 elements after the given sequence.
    """
    
    original_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]
    
    # We will generate the full sequence programmatically
    full_sequence = []
    
    # Initial state for the value generator
    value = 3
    direction = -1
    
    # Run index to determine the count
    run_index = 1
    
    # Generate enough elements to find the solution
    while len(full_sequence) < len(original_sequence) + 4:
        # Determine the count for the current run
        count = 1 if run_index <= 4 else 3
        
        # Add the run of numbers to the sequence
        for _ in range(count):
            full_sequence.append(value)
            
        # Update the value for the next run (bounce logic)
        if value == 1:
            direction = 1
        elif value == 3:
            direction = -1
        value += direction
        
        run_index += 1
        
    # Extract the next 4 elements
    next_4_elements = full_sequence[len(original_sequence):len(original_sequence) + 4]
    
    # The prompt asks to output each number in the "final equation".
    # We interpret this as showing how the original sequence is completed.
    completed_sequence = original_sequence + next_4_elements
    
    print("Original sequence: ", " ".join(map(str, original_sequence)))
    print("Next 4 elements are:", " ".join(map(str, next_4_elements)))
    print("Completed sequence:", " ".join(map(str, completed_sequence)))

solve_sequence()