import sys

def solve_ballet_query():
    """
    Analyzes classical ballet steps to find which one has the same starting and ending leg position.
    """
    # The definitions of the ballet steps and their properties are encoded here.
    # True: Ending position is the same as the starting position.
    # False: Ending position is different from the starting position.
    ballet_steps = {
        'A': {'name': 'Entrechat six', 'returns_to_start': True},
        'B': {'name': 'Échappé battu changé', 'returns_to_start': False},
        'C': {'name': 'Assemblé', 'returns_to_start': False},
        'D': {'name': 'Glissade derrière', 'returns_to_start': False},
        'E': {'name': 'Gargouillade', 'returns_to_start': False}
    }

    print("Analyzing ballet steps based on their starting and ending positions...")
    
    found_step = None
    for key, properties in ballet_steps.items():
        if properties['returns_to_start']:
            found_step = properties['name']
            break
            
    if found_step:
        print(f"\nThe step that ends in the same leg position as it started is: {found_step}")
        print("\nExplanation:")
        print("An Entrechat six is a vertical jump where the legs cross three times (making six movements). Even-numbered entrechats (like quatre, six, huit) always land in the same fifth position they started from.")
    else:
        print("Could not identify a step from the list that meets the criteria.")

solve_ballet_query()