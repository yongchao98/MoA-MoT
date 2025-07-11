def find_definitive_species():
    """
    Identifies the species in the image that definitively indicates an Atlantic,
    rather than Pacific, coast location.

    The image shows a community of organisms typical of the North Atlantic rocky shore.
    This includes:
    - Barnacles, identified as Semibalanus balanoides.
    - Snails, including the common periwinkle, Littorina littorea.

    While Littorina littorea is overwhelmingly Atlantic, it does have a very minor
    introduced presence on the Pacific coast.

    In contrast, the northern rock barnacle, Semibalanus balanoides, is found
    exclusively in the North Atlantic Ocean. Its presence is therefore a definitive
    indicator that the location is on an Atlantic coast.
    """
    definitive_species_name = "Semibalanus balanoides"
    print(definitive_species_name)

find_definitive_species()