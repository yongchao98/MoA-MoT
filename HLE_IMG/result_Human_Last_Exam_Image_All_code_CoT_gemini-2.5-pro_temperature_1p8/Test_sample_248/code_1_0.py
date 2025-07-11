def solve():
    """
    This function identifies the species from the image that is exclusively found on the Atlantic Coast
    and prints its scientific name.
    """
    # The image contains several species, including barnacles and two types of periwinkle snails.
    # The larger snail is the Common Periwinkle, Littorina littorea. While native to the Atlantic,
    # it has been introduced to the Pacific coast, so it's not a definitive indicator.
    # The smaller, darker snails are Rough Periwinkles, Littorina saxatilis.
    # The geographic range of Littorina saxatilis is restricted to the North Atlantic Ocean.
    # Therefore, its presence definitively indicates the photo was taken on an Atlantic coast.
    scientific_name = "Littorina saxatilis"
    print(scientific_name)

solve()