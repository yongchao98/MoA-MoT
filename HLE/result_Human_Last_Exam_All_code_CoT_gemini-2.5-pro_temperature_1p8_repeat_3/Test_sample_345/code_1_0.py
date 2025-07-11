import sys
import json
try:
    from SPARQLWrapper import SPARQLWrapper, JSON
except ImportError:
    print("Please install SPARQLWrapper: pip install SPARQLWrapper")
    sys.exit(1)

def find_museum_for_painting():
    """
    This function queries Wikidata to find the museum that acquired
    Kurt Günther's 'The Radionist' in 1967.
    """
    endpoint_url = "https://query.wikidata.org/sparql"
    # SPARQL query to find the painting, filter by collection and acquisition date.
    # We search for "Der Radionist", the original German title.
    # wdt:P170 is the property for "creator". wd:Q1794270 is the Wikidata item for "Kurt Günther".
    # rdfs:label is the name/title.
    # p:P195 is the property for "collection".
    # pq:P580 is the property for "start time" which is used for acquisition date here.
    query = """
    SELECT ?museumLabel WHERE {
      # Find the painting by its German title and creator
      ?painting rdfs:label "Der Radionist"@de;
                wdt:P170 wd:Q1794270.

      # Get the statement about the collection
      ?painting p:P195 ?collectionStatement.
      
      # Get the museum from the statement
      ?collectionStatement ps:P195 ?museum.
      
      # Get the acquisition date from the statement's qualifiers
      ?collectionStatement pq:P580 ?acquisitionDate.

      # Filter for the specific year 1967
      FILTER(YEAR(?acquisitionDate) = 1967)

      # Get the label (name) of the museum in English or German
      SERVICE wikibase:label { bd:serviceParam wikibase:language "en,de". }
    }
    LIMIT 1
    """

    sparql = SPARQLWrapper(endpoint_url)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)

    try:
        results = sparql.query().convert()
        bindings = results["results"]["bindings"]

        if bindings:
            museum_name = bindings[0]["museumLabel"]["value"]
            year_painted = 1927
            year_acquired = 1967
            painting_name = "The Radionist"
            
            # Print out the final answer with all the relevant numbers
            print(f"The {year_painted} tempera painting \"{painting_name}\" was acquired in the year {year_acquired} by the following museum:")
            print(museum_name)
        else:
            print("Could not find the museum based on the query criteria.")

    except Exception as e:
        print(f"An error occurred during the query: {e}")

if __name__ == "__main__":
    find_museum_for_painting()