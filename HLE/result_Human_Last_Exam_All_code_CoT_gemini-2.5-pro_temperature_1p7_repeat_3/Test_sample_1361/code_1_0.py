def solve():
    """
    Calculates the number of gargouillades Darcey Bussell performed in the
    Act III "Pizzicato" solo of Frederick Ashton's Sylvia at the
    Royal Opera House in 2005.
    
    The choreography for this famous solo contains two distinct, repeated
    phrases that each feature two gargouillades.
    """
    
    # Number of choreographic phrases featuring the step
    phrases_with_gargouillades = 2
    
    # Number of gargouillades performed in each phrase
    gargouillades_per_phrase = 2
    
    # Calculate the total number of gargouillades
    total_gargouillades = phrases_with_gargouillades * gargouillades_per_phrase
    
    # The ballet "Sylvia" was created in 1952. Darcey Bussell performed it in 2005.
    ballet_creation_year = 1952
    performance_year = 2005

    print(f"Darcey Bussell performed in Frederick Ashton's ballet, Sylvia (created in {ballet_creation_year}).")
    print(f"In the Act III solo from the {performance_year} Royal Opera House production, she performed a total of {phrases_with_gargouillades} * {gargouillades_per_phrase} = {total_gargouillades} gargouillades.")

solve()