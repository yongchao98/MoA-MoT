import numpy as np

def analyze_flag_rank(name, matrix_list):
    """
    Calculates and prints the rank of a flag's matrix representation.
    The final equation for the rank is printed by showing the matrix.
    """
    matrix = np.array(matrix_list)
    rank = np.linalg.matrix_rank(matrix)
    
    print(f"--- Analyzing: {name} ---")
    print("Matrix Representation:")
    for row in matrix_list:
        # Each number in the final equation is printed here
        print(f"  {row}")
    print(f"Result: The rank for {name} is {rank}.\n")
    return rank

def find_flags_with_matching_rank():
    """
    Models flags as matrices, calculates their ranks, and finds matches to Denmark's rank.
    To maximize rank, we assign distinct numerical values to distinct colors or symbols.
    Example mapping: Red=1, White=2, Blue=3, Green=4, Yellow=5, Black=6, Symbol/Other=7
    """
    
    # 1. Analyze Denmark to establish the target rank.
    # The Danish flag has two row patterns: [Red, White, Red] and [White, White, White].
    denmark_matrix = [
        [1, 2, 1],
        [2, 2, 2],
        [1, 2, 1]
    ]
    target_rank = analyze_flag_rank("Denmark", denmark_matrix)
    print(f"The target rank, based on the flag of Denmark, is {target_rank}.\n")
    
    # 2. Define matrix representations for various African flags.
    african_flags = {
        "Algeria": {
            "description": "A vertical bicolor (Green/White) with a central Red symbol. A row through the symbol is distinct from a row above or below it.",
            "matrix": [
                [4, 4, 2, 2],  # Row above symbol: [Green, Green, White, White]
                [4, 7, 7, 2],  # Row through symbol: [Green, Symbol, Symbol, White]
                [4, 4, 2, 2]   # Row below symbol
            ]
        },
        "Burundi": {
            "description": "A white saltire (diagonal cross) dividing fields of Red and Green, with a central circle symbol.",
            "matrix": [
                [4, 2, 1],  # Top section: [Green, White, Red]
                [2, 7, 2],  # Middle (cross/circle): [White, Symbol, White]
                [1, 2, 4]   # Bottom section: [Red, White, Green]
            ]
        },
        "Liberia": {
            "description": "A canton (Blue square) on a field of Red and White stripes. Rows within the canton are distinct from rows below it.",
            "matrix": [
                [3, 1],      # Row in canton (Blue), on a Red stripe
                [3, 2],      # Row in canton (Blue), on a White stripe
                [1, 1],      # Row below canton, on a Red stripe
                [2, 2]       # Row below canton, on a White stripe
            ]
        },
        "Togo": {
            "description": "A canton (Red square) on a field of Green and Yellow stripes. This has a similar structure to Liberia.",
            "matrix": [
                [1, 4],      # Row in canton (Red), on a Green stripe
                [1, 5],      # Row in canton (Red), on a Yellow stripe
                [4, 4],      # Row below canton, on a Green stripe
                [5, 5]       # Row below canton, on a Yellow stripe
            ]
        },
        "Nigeria": {
            "description": "A simple vertical tricolor (Green, White, Green). All rows are identical.",
            "matrix": [
                [4, 2, 4],
                [4, 2, 4],
                [4, 2, 4]
            ]
        },
        "Botswana": {
            "description": "Simple horizontal stripes (Blue, White, Black, White, Blue). All columns are identical, and rows are scalar multiples.",
            "matrix": [
                [3, 3, 3],  # Blue
                [2, 2, 2],  # White
                [6, 6, 6]   # Black
            ]
        },
        "South Africa": {
            "description": "A complex 'Y' shape design with 6 colors, resulting in at least 3 distinct row patterns.",
            "matrix": [
                [1, 2, 6],  # Top (Red, White, Black)
                [4, 5, 4],  # Middle (Green, Yellow, Green)
                [3, 2, 6]   # Bottom (Blue, White, Black)
            ]
        }
    }
    
    print("--- Analyzing African Flags ---")
    matching_flags = []
    for name, data in african_flags.items():
        print(f"Description: {data['description']}")
        rank = analyze_flag_rank(name, data['matrix'])
        if rank == target_rank:
            matching_flags.append(name)
            
    print("=========================================================")
    print(f"The African nations whose flags have the same rank ({target_rank}) as Denmark are:")
    for name in matching_flags:
        print(f"- {name}")
    print("=========================================================")

if __name__ == '__main__':
    find_flags_with_matching_rank()