def solve():
    """
    Identifies the animal in the image and provides its scientific name.
    """
    animal_present = True
    common_name = "Alligator Gar"
    scientific_name = "Atractosteus spatula"

    if animal_present:
        print(f"Yes, there is an animal in the image. It appears to be an {common_name}.")
        print(f"Its scientific name is {scientific_name}.")
    else:
        print("No animal could be clearly identified in the image.")

solve()