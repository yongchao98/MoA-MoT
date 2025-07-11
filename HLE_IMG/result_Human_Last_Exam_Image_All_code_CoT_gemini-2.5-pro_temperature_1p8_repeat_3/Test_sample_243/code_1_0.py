import wikipedia

def find_endemic_island():
    """
    Identifies the plant and finds the island it is endemic to.
    The plant is Procris boninensis. The name "boninensis" indicates its origin.
    """
    plant_name = "Procris boninensis"
    # The epithet 'boninensis' refers to the Bonin Islands.
    island_name = "Bonin Islands"

    print(f"The plant identified is likely {plant_name}.")
    print(f"This plant is endemic to the {island_name}.")
    
    # Let's use a library to confirm details about the Bonin Islands.
    try:
        page = wikipedia.page("Bonin Islands", auto_suggest=False)
        print("\nFrom Wikipedia:")
        print(page.summary.split('\n')[0]) # Print the first paragraph of the summary
    except wikipedia.exceptions.PageError:
        print("\nCould not retrieve additional information for the Bonin Islands from Wikipedia.")
    except Exception as e:
        # Fallback if Wikipedia is unreachable
        print("\nThe Bonin Islands (Ogasawara Islands) are an archipelago belonging to Japan.")


find_endemic_island()
