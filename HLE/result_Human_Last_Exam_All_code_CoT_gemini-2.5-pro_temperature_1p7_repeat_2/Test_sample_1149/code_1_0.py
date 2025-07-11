import sys

def get_anatomical_data(species, feature):
    """
    A conceptual function that retrieves specific anatomical data from a paleontological database.
    For this query, the information is based on the Sues et al. (2002) description of Irritator challengeri.
    """
    database = {
        "Irritator challengeri": {
            "lacrimal recess depressions": 2
        }
    }
    
    try:
        return database[species][feature]
    except KeyError:
        return "Data not found"

# Define the species and the anatomical feature we are looking for
species_name = "Irritator challengeri"
feature_name = "lacrimal recess depressions"

# Retrieve the number of depressions
number_of_depressions = get_anatomical_data(species_name, feature_name)

if isinstance(number_of_depressions, int):
    # Print the result in a full sentence as requested
    print(f"The number of smaller depressions contained within the larger lacrimal recess in {species_name} is: {number_of_depressions}")
else:
    print(number_of_depressions, file=sys.stderr)
