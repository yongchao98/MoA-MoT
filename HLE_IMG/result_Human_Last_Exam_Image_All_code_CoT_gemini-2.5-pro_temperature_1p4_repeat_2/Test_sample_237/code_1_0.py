def identify_animal():
    """
    Identifies the animal in the image and provides its scientific name.
    """
    is_animal_present = True
    common_name = "Longhorn Beetle"
    # Based on the morphology (long antennae, body shape, and coloration),
    # the insect is likely a male Stenocorus meridianus.
    scientific_name = "Stenocorus meridianus"

    if is_animal_present:
        print("Yes, there is an animal in this image.")
        print(f"It appears to be a type of {common_name}.")
        print(f"The likely scientific name is: {scientific_name}")
    else:
        print("No animal could be identified in the image.")

identify_animal()