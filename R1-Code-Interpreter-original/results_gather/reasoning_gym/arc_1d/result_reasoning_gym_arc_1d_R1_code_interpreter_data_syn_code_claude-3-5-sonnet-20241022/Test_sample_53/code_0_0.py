def analyze_shifts(input_str, output_str):
    input_list = [int(x) for x in input_str.split()]
    output_list = [int(x) for x in output_str.split()]
    
    # Find differences
    print("Length:", len(input_list))
    print("Input:", input_list)
    print("Output:", output_list)
    
    # Find where values differ
    differences = [(i, input_list[i], output_list[i]) 
                  for i in range(len(input_list)) 
                  if input_list[i] != output_list[i]]
    print("Differences at positions:", differences)

# Example 1
print("Example 1:")
analyze_shifts("6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 0 0 6",
              "6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 0 0 6 6 6")

print("\nExample 3:")
analyze_shifts("6 6 6 6 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 6 6 6 6 6 6 6 6 6",
              "6 6 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 6 6 6 6 6 6 6 6 6 6 6")