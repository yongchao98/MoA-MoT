def analyze_movement(input_str, output_str):
    input_grid = [int(x) for x in input_str.split()]
    output_grid = [int(x) for x in output_str.split()]
    
    # Find positions of non-zero numbers
    def find_blocks(grid):
        blocks = []
        current_block = []
        current_num = 0
        
        for i, num in enumerate(grid):
            if num != 0:
                if num == current_num:
                    current_block.append(i)
                else:
                    if current_block:
                        blocks.append((current_num, current_block))
                    current_block = [i]
                    current_num = num
            else:
                if current_block:
                    blocks.append((current_num, current_block))
                    current_block = []
                    current_num = 0
        if current_block:
            blocks.append((current_num, current_block))
        return blocks

    print("Example analysis:")
    print("Input blocks:", find_blocks(input_grid))
    print("Output blocks:", find_blocks(output_grid))

# Analyze all examples
print("Example 1:")
analyze_movement("0 0 0 0 0 0 0 0 0 2 0 6 6 6 0 0", "0 0 0 0 0 6 6 6 0 2 0 0 0 0 0 0")
print("\nExample 2:")
analyze_movement("0 0 0 0 0 4 4 2 0 0 0 0 0 0 0 0", "0 0 0 0 0 0 0 2 4 4 0 0 0 0 0 0")
print("\nExample 3:")
analyze_movement("0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 9", "0 0 0 9 0 0 0 0 0 2 0 0 0 0 0 0")