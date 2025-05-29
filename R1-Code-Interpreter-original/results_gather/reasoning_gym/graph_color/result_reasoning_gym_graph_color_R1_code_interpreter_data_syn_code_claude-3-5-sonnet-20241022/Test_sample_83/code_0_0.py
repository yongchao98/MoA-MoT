import json

def color_graph():
    # Initialize all vertices with color 1
    coloring = {str(i): 1 for i in range(10)}
    
    # Handle the only edge (1,5)
    coloring['1'] = 1
    coloring['5'] = 2
    
    # Convert to JSON string
    result = json.dumps(coloring)
    print(result)

color_graph()