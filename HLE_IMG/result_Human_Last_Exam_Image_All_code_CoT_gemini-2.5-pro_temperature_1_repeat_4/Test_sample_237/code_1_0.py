def identify_animal():
    """
    Identifies the animal in the image and provides its scientific name.
    """
    is_animal_present = True
    animal_type = "Longhorn Beetle"
    scientific_name = "Stictoleptura cordigera"

    if is_animal_present:
        print("Yes, there is an animal in this image.")
        print(f"It appears to be a {animal_type}.")
        print(f"The animal's scientific name is: {scientific_name}")

identify_animal()