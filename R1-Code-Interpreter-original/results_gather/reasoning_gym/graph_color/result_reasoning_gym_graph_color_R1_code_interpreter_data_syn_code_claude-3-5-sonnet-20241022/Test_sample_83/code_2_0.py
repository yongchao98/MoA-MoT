import json

def color_graph():
    # Initialize all vertices with color 1
    coloring = {str(i): 1 for i in range(10)}
    
    # For the edge (1,5):
    # Let's explicitly make them different colors
    coloring['1'] = 1
    coloring['5'] = 2  # Changed this to ensure different color from vertex 1
    
    # Convert to JSON string
    result = json.dumps(coloring)
    print(result)

color_graph()