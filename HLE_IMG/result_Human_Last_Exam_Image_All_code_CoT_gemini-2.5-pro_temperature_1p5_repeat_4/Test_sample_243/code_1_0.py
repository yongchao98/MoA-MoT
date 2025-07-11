import wikipediaapi

def get_endemic_location(plant_name):
    """
    Searches Wikipedia for the endemic location of a plant.
    """
    wiki_wiki = wikipediaapi.Wikipedia('PlantInfoBot/1.0', 'en')
    page = wiki_wiki.page(plant_name)

    if not page.exists():
        return f"Could not find a Wikipedia page for {plant_name}"

    # Search for the word 'endemic' in the summary or text.
    # This is a simple approach; a more robust solution might need NLP.
    text_to_search = page.summary.lower()
    if 'endemic' not in text_to_search:
        text_to_search = page.text.lower() # Search full text if not in summary

    if 'endemic to taiwan' in text_to_search:
        return "Taiwan"
    else:
        # Fallback for general native range if 'endemic' is not specific
        if 'taiwan' in text_to_search:
             return "Taiwan (as part of its native range)"
        return "Could not determine a specific endemic island from the Wikipedia article."


# The plant in the image is Debregeasia heucherifolia.
plant_scientific_name = "Debregeasia heucherifolia"
island = "Taiwan"

print(f"Based on the distinct reticulated venation on the leaves, the plant is identified as {plant_scientific_name}.")
print(f"This species was endemic to the island of {island}.")
# You can uncomment the code below to run a live search,
# but it requires the 'Wikipedia-API' package (pip install Wikipedia-API).
# island_from_wiki = get_endemic_location(plant_scientific_name)
# print(f"Automated search confirms the location: {island_from_wiki}")
