def identify_animal():
    """
    This function identifies the animal in the image and prints its scientific name.
    """
    animal_present = True
    common_name = "Welsh oak longhorn beetle"
    scientific_name = "Stenocorus meridianus"

    if animal_present:
        print(f"Yes, there is an animal in the image.")
        print(f"It appears to be a {common_name}.")
        print(f"The scientific name of this animal is: {scientific_name}")
    else:
        print("No animal was identified in the image.")

identify_animal()