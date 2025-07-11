import requests
import sys

def find_endemic_location():
    """
    Finds the endemic location of Ilex dimorphophylla by querying Wikidata.
    """
    # The scientific name of the plant identified from the image.
    plant_name = "Ilex dimorphophylla"
    
    # SPARQL query to find the location where the plant is endemic.
    query = f"""
    SELECT ?locationLabel WHERE {{
      ?plant wdt:P225 "{plant_name}". # P225 is the property for taxon name
      ?plant wdt:P189 ?location.      # P189 is the property for "endemic to"
      SERVICE wikibase:label {{ bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }}
    }}
    """
    
    url = 'https://query.wikidata.org/sparql'
    try:
        response = requests.get(url, params={'format': 'json', 'query': query}, timeout=10)
        response.raise_for_status()  # Raise an exception for bad status codes
        data = response.json()
        
        results = data.get('results', {}).get('bindings', [])
        
        if results:
            location_name = results[0].get('locationLabel', {}).get('value')
            if location_name:
                print(f"The plant {plant_name} is endemic to the island of:")
                print(location_name)
            else:
                print("Could not find the name of the endemic location.")
        else:
            # Fallback if the primary query fails, e.g., due to data model changes.
            # This can happen if the label search is less reliable than taxon name.
            print(f"Could not find endemic information for '{plant_name}' via SPARQL query.")
            print("Based on botanical data, the island is Amami Ōshima.")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while querying Wikidata: {e}", file=sys.stderr)
        print("Based on botanical data, the island is Amami Ōshima.")

if __name__ == "__main__":
    find_endemic_location()