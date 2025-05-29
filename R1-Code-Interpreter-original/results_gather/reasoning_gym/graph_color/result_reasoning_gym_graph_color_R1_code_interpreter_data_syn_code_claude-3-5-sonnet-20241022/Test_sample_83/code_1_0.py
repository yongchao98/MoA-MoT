import json

def color_graph():
    # Initialize all vertices with color 1
    coloring = {str(i): 1 for i in range(10)}
    
    # Handle the only edge (1,5)
    # Vertex 1 gets color 2, Vertex 5 gets color 1
    coloring['1'] = 2  # Changed this from 1 to 2
    coloring['5'] = 1
    
    # Convert to JSON string
    result = json.dumps(coloring)
    print(result)

color_graph()