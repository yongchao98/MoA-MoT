import json

def verify_coloring():
    # Initialize vertices and edges
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    edges = [(1, 8), (7, 8)]
    
    # Create the coloring
    coloring = {
        0: 1,  # isolated vertex
        1: 2,  # connected to 8
        2: 1,  # isolated vertex
        3: 1,  # isolated vertex
        4: 1,  # isolated vertex
        5: 1,  # isolated vertex
        6: 1,  # isolated vertex
        7: 2,  # connected to 8
        8: 1,  # connected to 1 and 7
        9: 1   # isolated vertex
    }
    
    # Verify the coloring
    for v1, v2 in edges:
        if coloring[v1] == coloring[v2]:
            print("Invalid coloring!")
            return None
    
    # Convert to JSON string format
    result = json.dumps(coloring)
    print(result)

verify_coloring()