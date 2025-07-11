import numpy as np

def solve():
    """
    Solves the grid transformation puzzle.
    """
    # The test input grid for the puzzle
    test_input = np.array([
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ])

    # The 3x3 kernel for the convolution, determined from the examples
    kernel = np.array([
        [1, 0, 0],
        [1, 1, 1],
        [0, 1, 0]
    ])
    
    size = len(test_input)
    output_grid = np.zeros((size, size), dtype=int)
    
    # Apply the kernel to each cell of the input grid
    for r in range(size):
        for c in range(size):
            total = 0
            # Iterate through the 3x3 neighborhood
            for dr in range(-1, 2):
                for dc in range(-1, 2):
                    # Kernel value for this position
                    k_val = kernel[dr + 1, dc + 1]
                    if k_val == 1:
                        # Get input grid value with wrap-around (toroidal)
                        in_val = test_input[(r + dr) % size, (c + dc) % size]
                        total += in_val
            
            # The output is the total sum modulo 2
            output_grid[r, c] = total % 2
            
    # As requested, print the numbers in the final equation for each output cell.
    # The equation is out(r,c) = (sum of neighbors) mod 2. We print the values of the neighbors.
    print("Applying transformation out(r,c) = (in(r-1,c-1) + in(r,c-1) + in(r,c) + in(r,c+1) + in(r+1,c)) % 2")
    for r in range(size):
        for c in range(size):
            # These are the neighbors included in the sum based on the kernel [1,0,0], [1,1,1], [0,1,0]
            val1 = test_input[(r - 1) % size, (c - 1) % size]
            val2 = test_input[(r) % size, (c - 1) % size]
            val3 = test_input[(r) % size, (c) % size]
            val4 = test_input[(r) % size, (c + 1) % size]
            val5 = test_input[(r + 1) % size, (c) % size]
            result = (val1 + val2 + val3 + val4 + val5) % 2
            print(f"out({r},{c}) = ({val1} + {val2} + {val3} + {val4} + {val5}) % 2 = {result}")

    # Flatten the grid and join to form the final binary string
    final_string = "".join(map(str, output_grid.flatten()))
    print("\nFinal Output Grid:")
    print(output_grid)
    print("\nFinal Answer String:")
    print(final_string)

solve()