def solve_trilogy():
    """
    This function solves the three-part riddle to find a trilogy of cities.
    """
    
    # 1. In Irving Wallace's "The Miracle", a young girl anticipates an encounter with a lady in white. 
    # This refers to the Marian apparitions to Bernadette Soubirous.
    # The city where this took place is Lourdes, France.
    city_1 = "Lourdes"
    
    # 2. A poem by Yuri Krasnokutsky describes a cold stream flowing through the blood of 'X'. 
    # The poem is about Lenin ("X"). Lenin's Mausoleum is in Red Square.
    # The city is Moscow.
    city_2 = "Moscow"
    
    # 3. The Croatian gastronomic paradise "Dolac" (in Zagreb) is nicknamed the "Belly of Zagreb". 
    # Its most famous analogue is "The Belly of Paris" (Le Ventre de Paris), a term popularized by
    # Émile Zola's novel about the Les Halles market.
    # The city is Paris.
    city_3 = "Paris"
    
    # Print the final trilogy in the specified format.
    print(f"{{{city_1}, {city_2}, {city_3}}}")

solve_trilogy()