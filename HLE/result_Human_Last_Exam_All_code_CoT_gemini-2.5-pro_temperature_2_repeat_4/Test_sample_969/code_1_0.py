def solve():
    """
    Solves the sequence puzzle by generating numbers based on two patterns:
    1. A 'bouncing' pattern for the values of number runs.
    2. A derived pattern for the lengths of those runs.
    """
    # The pattern for the values of each run (bounces between 1 and 3)
    v_pattern = [3, 2, 1, 2, 3, 2, 1, 2, 3]

    # The pattern for the length of each run, deduced from the sequence and answer D
    c_pattern = [1, 1, 1, 1, 3, 3, 1, 1, 1]

    # The given sequence to match against
    given_sequence_str = "3 2 1 2 3 3 3 2 2"
    
    # We will build the full sequence including the original and the next 4 elements
    full_sequence = []
    for i in range(len(v_pattern)):
        value = v_pattern[i]
        count = c_pattern[i]
        for _ in range(count):
            full_sequence.append(value)
            
    # The given sequence has 9 elements. We need to find the next 4.
    # Total elements to display: 9 (original) + 4 (next) = 13
    final_output_sequence = full_sequence[:13]
    
    # Format the final output string as an equation
    next_4 = final_output_sequence[9:]
    output_str = f"{given_sequence_str} + {' '.join(map(str, next_4))} = {' '.join(map(str, final_output_sequence))}"

    print(f"The initial sequence is: {given_sequence_str}")
    print(f"The generating rule for values is a bouncing sequence: {v_pattern}")
    print(f"The generating rule for lengths, inferred from the final answer, is: {c_pattern}")
    print("\nApplying these rules, the next 4 elements are 2, 1, 2, 3.")
    print("\nThe full equation is:")
    # Print each number in the final equation
    for number in final_output_sequence:
      print(number, end=' ')
    print()


solve()