def identify_animal():
    """
    Identifies the animal in the image and provides its scientific name.
    """
    animal_present = True
    animal_name = "Moose"
    scientific_name = "Alces alces"
    body_part = "snout (muffle)"

    if animal_present:
        print(f"Yes, there is an animal in the image. It appears to be the {body_part} of a {animal_name}.")
        print(f"The scientific name of this animal is: {scientific_name}")
    else:
        print("No animal could be clearly identified in the image.")

identify_animal()