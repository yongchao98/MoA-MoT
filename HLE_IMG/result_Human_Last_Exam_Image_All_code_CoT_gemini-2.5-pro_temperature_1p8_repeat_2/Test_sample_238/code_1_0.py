def identify_organism():
    """
    Identifies the small black organisms in the image.
    These are the tiny, hair-like structures with black heads growing on top of the
    larger, white and brown shelf fungi.
    """
    
    # The organism is identified as a slime mold, specifically its reproductive stage.
    organism_name = "Slime mold (Myxomycete)"
    common_example = "Comatricha nigra"
    description = "fruiting bodies"
    
    print(f"The small, black organisms visible on the surface of the larger fungus are the {description} of a {organism_name}.")
    print(f"These are not true fungi but amoeba-like protists. A common species with this appearance is {common_example}, often called hairy stereum slime.")

identify_organism()