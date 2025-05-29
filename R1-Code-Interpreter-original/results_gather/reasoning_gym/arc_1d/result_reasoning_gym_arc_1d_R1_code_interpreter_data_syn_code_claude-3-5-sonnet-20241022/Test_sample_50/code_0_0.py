def analyze_pattern(input_str, output_str):
    input_arr = [int(x) for x in input_str.split()]
    output_arr = [int(x) for x in output_str.split()]
    
    # Find non-zero number
    num = next(x for x in input_arr if x != 0)
    
    # Count numbers and positions
    input_count = sum(1 for x in input_arr if x == num)
    output_count = sum(1 for x in output_arr if x == num)
    
    # Find first and last positions
    input_first = input_arr.index(num)
    input_last = len(input_arr) - 1 - input_arr[::-1].index(num)
    output_first = output_arr.index(num)
    output_last = len(output_arr) - 1 - output_arr[::-1].index(num)
    
    print(f"Number: {num}")
    print(f"Count: Input={input_count}, Output={output_count}")
    print(f"First position: Input={input_first}, Output={output_first}")
    print(f"Last position: Input={input_last}, Output={output_last}")
    print("---")

# Analyze each example
print("Example 1:")
analyze_pattern("3 0 0 0 0 0 0 0 0 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3", 
                "3 3 3 3 0 0 0 0 0 0 0 0 3 3 3 3 3 3 3 3 3 3 3 3 3")

print("Example 2:")
analyze_pattern("5 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 5 5",
                "5 5 5 5 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5")

print("Example 3:")
analyze_pattern("1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1",
                "0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0")