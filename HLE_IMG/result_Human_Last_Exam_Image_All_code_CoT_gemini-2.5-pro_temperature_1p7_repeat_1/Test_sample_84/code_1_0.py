def find_city_with_similar_coat_of_arms():
    """
    This function identifies and prints the name of the city with a coat of arms
    that is compositionally very similar to the one in the provided image.
    
    The analysis is as follows:
    - The image shows a shield divided horizontally (per fess).
    - The top part is blue (azure) and contains a crescent.
    - The bottom part is gold (or) and contains a crab.
    
    A search for real-world coats of arms reveals that the municipal coat of arms
    of Ílhavo, Portugal, has an identical structure.
    
    - Ílhavo's shield: Divided horizontally, blue over gold.
    - Ílhavo's top part: Contains a crescent (officially silver, not gold).
    - Ílhavo's bottom part: Contains a crab (officially red, not black).
    
    The fundamental design is the same, with the differences in color and external
    ornaments (the crown and branches) being artistic variations.
    Therefore, Ílhavo is the city with the most similar coat of arms.
    """
    city_name = "Ílhavo"
    print(city_name)

find_city_with_similar_coat_of_arms()