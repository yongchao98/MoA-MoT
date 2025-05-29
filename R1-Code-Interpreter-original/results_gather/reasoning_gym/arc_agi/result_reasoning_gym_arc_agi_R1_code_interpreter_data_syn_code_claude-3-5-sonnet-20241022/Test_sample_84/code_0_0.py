def analyze_examples():
    # Example inputs and outputs
    examples = [
        {
            'input': [
                [9,9,9,9,6,9,9],
                [9,9,6,6,6,9,9],
                [6,6,6,9,9,9,9],
                [9,9,9,9,9,9,9],
                [9,9,9,6,6,6,9],
                [6,6,6,6,9,6,6],
                [9,9,9,9,9,9,9]
            ],
            'output': [9,9,9]
        },
        {
            'input': [
                [9,9,6,9,9,9,9,9,9],
                [9,9,6,9,9,9,9,9,9],
                [6,6,6,9,9,9,9,9,9],
                [9,9,9,9,9,9,6,6,6],
                [9,9,9,9,6,6,6,9,9],
                [6,6,6,6,6,9,9,9,9],
                [9,9,9,9,9,9,6,6,6],
                [9,6,6,6,6,6,6,9,9],
                [6,6,9,9,9,9,9,9,9],
                [6,9,9,6,6,6,6,9,9],
                [6,9,9,6,9,9,6,6,6],
                [9,9,9,6,9,9,9,9,9]
            ],
            'output': [9,9,9,9,9]
        },
        {
            'input': [
                [6,6,9],
                [9,6,6],
                [9,9,9]
            ],
            'output': [9,9]
        }
    ]
    
    for idx, example in enumerate(examples):
        input_grid = example['input']
        output_length = len(example['output'])
        rows = len(input_grid)
        cols = len(input_grid[0])
        
        # Count number of connected regions of 6's
        visited = [[False] * cols for _ in range(rows)]
        regions = 0
        
        def dfs(i, j):
            if i < 0 or i >= rows or j < 0 or j >= cols or visited[i][j] or input_grid[i][j] != 6:
                return
            visited[i][j] = True
            # Check all 8 directions
            for di in [-1,0,1]:
                for dj in [-1,0,1]:
                    dfs(i+di, j+dj)
        
        for i in range(rows):
            for j in range(cols):
                if not visited[i][j] and input_grid[i][j] == 6:
                    regions += 1
                    dfs(i, j)
        
        print(f"Example {idx+1}:")
        print(f"Grid size: {rows}x{cols}")
        print(f"Number of connected 6's regions: {regions}")
        print(f"Output length: {output_length}")
        print()

analyze_examples()