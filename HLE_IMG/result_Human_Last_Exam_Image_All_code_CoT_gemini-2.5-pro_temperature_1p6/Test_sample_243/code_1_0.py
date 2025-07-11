import wikipediaapi

def get_endemic_island_from_plant(plant_species_epithet):
    """
    Looks up the island associated with a plant's specific epithet.
    In this case, the logic is pre-determined from visual identification.
    """
    # The plant in the image is identified as a Helwingia species.
    # The species Helwingia formosana is endemic to a specific island.
    # The epithet 'formosana' refers to Formosa.
    if plant_species_epithet == "formosana":
        # Formosa is the historical name for Taiwan.
        island_name = "Taiwan"
        print(f"The plant Helwingia formosana is endemic to the island of {island_name}.")
    else:
        print("Could not determine the endemic island for the given plant.")
        island_name = None
    
    return island_name

# The plant's specific epithet suggesting an island origin is "formosana".
final_answer = get_endemic_island_from_plant("formosana")

# Although the function prints the sentence, the final requested output is just the name.
# print(final_answer) # This would be the most direct way to output.
# Since the prompt asks to output the final answer directly in the response,
# and not to ask the user to copy/paste, I will format it as requested.
