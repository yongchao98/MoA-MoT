def identify_animal():
    """
    Identifies the animal in the image and prints its scientific name.
    """
    animal_present = True
    common_name = "Pinecone fish"
    scientific_name = "Monocentris japonica"

    if animal_present:
        print(f"Yes, there is an animal in the image. It is a {common_name}.")
        print(f"Its scientific name is: {scientific_name}")
    else:
        print("No animal could be identified in the image.")

identify_animal()