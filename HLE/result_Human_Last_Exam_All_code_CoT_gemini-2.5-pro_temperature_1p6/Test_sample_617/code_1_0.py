def get_pseudomonas_color():
    """
    Determines the color of a concentrated Pseudomonas aeruginosa sample
    based on its known pigment production.
    """
    bacterium = "Pseudomonas aeruginosa"
    
    # Define the primary pigments and their colors
    pigment_1 = {"name": "Pyocyanin", "color": "Blue"}
    pigment_2 = {"name": "Pyoverdin", "color": "Green"}
    
    # In dense cultures, both pigments are produced and mix together.
    # The combination of Blue and Green gives a Blue-green color.
    resulting_color = "Blue-green"
    
    print(f"The organism is {bacterium}.")
    print(f"It produces the pigment {pigment_1['name']}, which has a {pigment_1['color']} color.")
    print(f"It also produces the pigment {pigment_2['name']}, which has a {pigment_2['color']} color.")
    print(f"When these pigments mix in a dense culture, the resulting color is: {resulting_color}.")
    print("\nWashing and concentrating the bacteria makes this color more intense.")
    print("\nTherefore, the correct answer is C.")

# Execute the function to display the answer
get_pseudomonas_color()