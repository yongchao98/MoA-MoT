def identify_organism_in_photo():
    """
    Identifies the organisms in the photo, focusing on the black hair-like structures.
    """
    
    # Information about the black organism inquired about
    black_organism = {
        "common_name": "Bonnet Mold",
        "scientific_name": "Spinellus fusiger",
        "description": "A type of parasitic pin mold that grows on the caps of other mushrooms."
    }

    # Information about the host fungus
    host_fungus = {
        "common_name": "Purple-toothed Polypore (probable)",
        "scientific_name": "Trichaptum biforme (or a related species)",
        "description": "The larger shelf fungus that the Bonnet Mold is growing on."
    }

    print("Identification of the black organism in the photo:")
    print("-------------------------------------------------")
    print(f"Common Name: {black_organism['common_name']}")
    print(f"Scientific Name: {black_organism['scientific_name']}")
    print(f"Description: {black_organism['description']}")
    print("\nThis mold is growing on a larger host fungus, which is likely:")
    print(f"Host Common Name: {host_fungus['common_name']}")
    print(f"Host Scientific Name: {host_fungus['scientific_name']}")

if __name__ == "__main__":
    identify_organism_in_photo()