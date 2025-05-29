import json

def color_graph():
    # Initialize all vertices with color 1
    coloring = {str(i): 1 for i in range(10)}
    
    # For the edge (1,5):
    # Vertex 1 gets color 2
    # Vertex 5 gets color 3
    # This ensures they have different colors
    coloring['1'] = 2
    coloring['5'] = 3
    
    # Convert to JSON string
    result = json.dumps(coloring)
    print(result)

color_graph()