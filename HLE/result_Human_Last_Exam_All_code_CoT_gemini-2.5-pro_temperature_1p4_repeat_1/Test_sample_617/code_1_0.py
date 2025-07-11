def get_sample_color():
    """
    Determines the color of a concentrated Pseudomonas aeruginosa sample.
    """
    organism = "Pseudomonas aeruginosa"
    
    # This bacterium is well-known for producing characteristic pigments.
    pigment1 = {
        "name": "Pyocyanin",
        "color": "Blue"
    }
    
    pigment2 = {
        "name": "Pyoverdine",
        "color": "Green" # or yellow-green
    }
    
    # A dense culture of P. aeruginosa will contain both of these pigments,
    # which are secreted by the cells. When the cells are concentrated into a
    # pellet, the combined color of the pigments becomes visible.
    
    combined_color = "Blue-green"
    
    print(f"The organism in question is {organism}.")
    print(f"It produces a pigment called {pigment1['name']}, which is {pigment1['color']}.")
    print(f"It also produces a pigment called {pigment2['name']}, which is {pigment2['color']}.")
    print("\nWhen these pigments are present together in a dense sample, they mix.")
    print(f"The resulting color of the concentrated sample is therefore: {combined_color}.")

# Run the function to display the answer.
get_sample_color()