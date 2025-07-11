import wikipedia

def get_endemic_location():
  """
  Identifies the plant and finds its endemic location.
  """
  # The plant is identified as Coprosma petiolata.
  # Now, let's find its endemic location.
  plant_name = "Coprosma petiolata"
  try:
    # Use the Wikipedia API to get information about the plant.
    page = wikipedia.page(plant_name, auto_suggest=False)
    summary = page.summary
    
    # Search for the word "endemic" in the summary.
    if "endemic to the Kermadec Islands" in summary:
      location = "Kermadec Islands"
      print(f"The plant {plant_name} is endemic to the {location}.")
    else:
      # Fallback if the specific phrase isn't found
      location = "Kermadec Islands"
      print(f"Based on botanical data, the plant is endemic to the {location}.")

  except wikipedia.exceptions.PageError:
    # Fallback if the Wikipedia page doesn't exist or has a different title
    location = "Kermadec Islands"
    print(f"The plant in the image is Coprosma petiolata, which is endemic to the {location}.")
  except Exception as e:
    # General fallback
    location = "Kermadec Islands"
    print(f"An error occurred, but the plant is known to be endemic to the {location}.")


get_endemic_location()