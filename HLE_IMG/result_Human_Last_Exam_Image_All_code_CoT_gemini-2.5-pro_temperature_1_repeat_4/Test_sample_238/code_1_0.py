def identify_organism():
    """
    This function identifies the small black organism seen growing on the
    larger fungus in the provided image.
    """
    # The main fungus is a type of bracket fungus (polypore).
    # The small black organisms growing on it are a parasitic mold.
    
    organism_info = {
        "common_name": "Bonnet Mould",
        "scientific_name": "Spinellus fusiger",
        "description": "A type of pin mold that is parasitic on the fruiting bodies of other fungi, especially those from the genus Mycena and Trametes."
    }

    print(f"The small, black, pin-like organism growing on the larger fungus is commonly known as {organism_info['common_name']}.")
    print(f"Its scientific name is {organism_info['scientific_name']}.")
    print(f"Description: {organism_info['description']}")

# Execute the function to print the identification.
identify_organism()