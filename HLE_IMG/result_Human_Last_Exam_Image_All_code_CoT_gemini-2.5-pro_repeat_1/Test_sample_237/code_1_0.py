def identify_animal():
    """
    Identifies the animal in the image and provides its scientific name.
    """
    is_animal_present = True
    common_name = "Red-brown Longhorn Beetle"
    scientific_name = "Stictoleptura rubra"

    if is_animal_present:
        print(f"Yes, there is an animal in the image.")
        print(f"It is a {common_name}.")
        print(f"The animal's scientific name is: {scientific_name}")
    else:
        print("There is no animal in the image.")

identify_animal()