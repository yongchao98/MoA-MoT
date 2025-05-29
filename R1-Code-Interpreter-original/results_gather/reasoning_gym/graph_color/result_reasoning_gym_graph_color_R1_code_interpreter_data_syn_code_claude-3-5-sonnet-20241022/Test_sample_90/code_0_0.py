import json

def get_graph_coloring():
    # Initialize all vertices with color 1 (any color would work for disconnected vertices)
    coloring = {str(i): 1 for i in range(10)}
    
    # Handle the connected components
    # Vertex 9 is connected to both 2 and 5
    # Let's give vertex 9 color 1
    coloring['9'] = 1
    
    # Then vertices 2 and 5 must have a different color than 9
    # They can have the same color since they're not connected
    coloring['2'] = 2
    coloring['5'] = 2
    
    # Convert to proper format and return
    result = json.dumps(coloring)
    print(result)

get_graph_coloring()