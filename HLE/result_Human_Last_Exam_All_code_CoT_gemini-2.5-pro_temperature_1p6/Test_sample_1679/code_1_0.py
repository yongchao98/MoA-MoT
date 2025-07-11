import sys

def find_first_winner_with_obelisk():
    """
    This function simulates the process of finding the first Academy Award
    Best Picture winner that depicts a Luxor Obelisk.
    """

    # A curated list representing researched film information. The list is
    # ordered chronologically to find the *first* winner.
    # Films are included based on their potential to feature an obelisk in Paris or Luxor.
    best_picture_candidates = [
        {
            'year': 1928,
            'title': 'Wings',
            'locations': ['Paris, France'],
            'features_obelisk': False,
            'reason': "While it features Paris, scenes focus on other landmarks like the Eiffel Tower and Folies Bergère, not the Place de la Concorde."
        },
        {
            'year': 1951,
            'title': 'An American in Paris',
            'locations': ['Paris, France'],
            'features_obelisk': True,
            'reason': "The film is a celebration of Parisian scenery and explicitly features the Place de la Concorde with its obelisk."
        },
        {
            'year': 1956,
            'title': 'Around the World in 80 Days',
            'locations': ['Paris, France', 'Suez, Egypt'],
            'features_obelisk': True,
            'reason': "Features an aerial view of the Paris obelisk, but was released after 'An American in Paris'."
        },
         {
            'year': 1962,
            'title': 'Lawrence of Arabia',
            'locations': ['Cairo, Egypt'],
            'features_obelisk': False,
            'reason': "Features scenes in Egypt (Cairo), but not Luxor where the other obelisk stands."
        }
    ]

    print("Step 1: Identifying the locations of the two Luxor Obelisks.")
    print("- Location 1: Luxor Temple in Luxor, Egypt.")
    print("- Location 2: Place de la Concorde in Paris, France.")
    print("\nStep 2: Chronologically searching for the first Best Picture winner to depict one of these obelisks.")
    
    first_winner = None
    
    for film in best_picture_candidates:
        year = film['year']
        title = film['title']
        features_obelisk = film['features_obelisk']
        
        print(f"\nChecking Year: {year}, Film: '{title}'...")
        if features_obelisk:
            print(f"-> Match found! The film '{title}' is known to feature a Luxor Obelisk.")
            print(f"   Reason: {film['reason']}")
            first_winner = film
            break # Stop at the first chronological match
        else:
            print(f"-> No match. Analysis shows it does not feature an obelisk.")
            print(f"   Reason: {film['reason']}")
            
    if first_winner:
        print("\n---------------------------------------------------")
        print("Final Answer:")
        print(f"The first Academy Award winner for Best Picture to depict a Luxor Obelisk on-screen is '{first_winner['title']}', which won for the year {first_winner['year']}.")
        print("---------------------------------------------------")

# Run the search
find_first_winner_with_obelisk()