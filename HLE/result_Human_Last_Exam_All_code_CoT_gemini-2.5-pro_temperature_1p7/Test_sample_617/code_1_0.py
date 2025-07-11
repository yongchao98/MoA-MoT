def get_pseudomonas_color():
    """
    Determines the color of a dense Pseudomonas aeruginosa sample
    based on its known pigment production.
    """
    # Pseudomonas aeruginosa is famous for producing several pigments.
    # The two most prominent ones that contribute to its visible color are:

    # 1. Pyocyanin, which has a blue color.
    pigment1 = {"name": "Pyocyanin", "color": "blue"}

    # 2. Pyoverdine (also known as fluorescein), which is a siderophore
    # with a yellow-green fluorescent color.
    pigment2 = {"name": "Pyoverdine", "color": "green"}

    print("Analysis of Pseudomonas aeruginosa sample color:")
    print(f"The bacterium produces {pigment1['name']}, which is {pigment1['color']}.")
    print(f"The bacterium also produces {pigment2['name']}, which is {pigment2['color']}.")
    
    # In a dense culture, both pigments are present. The washing and concentration
    # process results in a very dense pellet of cells, which will be strongly colored
    # by the combination of these pigments.
    combined_color = f"{pigment1['color']}-{pigment2['color']}"

    print(f"\nThe combination of '{pigment1['color']}' and '{pigment2['color']}' gives the sample a characteristic '{combined_color}' appearance.")
    print("\nTherefore, the dense, concentrated sample is blue-green.")

get_pseudomonas_color()