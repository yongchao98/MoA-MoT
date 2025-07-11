def identify_atlantic_species():
    """
    Identifies the species in the image that is characteristic of the
    Atlantic coast and not the Pacific coast.
    """
    # The image displays several intertidal organisms, including barnacles and snails.
    # The larger snail with spiral bands is the Common Periwinkle.
    # Scientific name: Littorina littorea
    # Distribution:
    # - Native to the Northeast Atlantic (Europe).
    # - Introduced and widespread on the Northwest Atlantic coast (North America).
    # - Not established on the Pacific coast of North America.
    # The other visible species, like the Northern Rock Barnacle (Semibalanus balanoides)
    # and the Rough Periwinkle (Littorina saxatilis), are found on both coasts.
    # Therefore, the presence of Littorina littorea is the key indicator.

    species_name = "Littorina littorea"
    
    print("The species that indicates this photo was taken on the Atlantic Coast is the Common Periwinkle.")
    print(f"Its scientific name is: {species_name}")

identify_atlantic_species()