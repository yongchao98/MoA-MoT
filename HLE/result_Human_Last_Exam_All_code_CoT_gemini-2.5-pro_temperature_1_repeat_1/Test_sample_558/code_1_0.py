def find_bass_singer():
    """
    This function identifies the bass singer from the described opera production.
    
    1. Opera Identification: Marietta Alboni's 1843 La Scala revival was Rossini's "Semiramide".
    2. Production Identification: The Metropolitan Opera's 1990 production of "Semiramide"
       matches the timeline (more than 70 years after Caruso's era, which ended in 1920).
    3. Role Identification: The principal bass role in "Semiramide" is Assur.
    4. Singer Identification: The bass who sang Assur in the Met's 1990 production was Samuel Ramey.
    """
    
    opera = "Semiramide"
    year_of_nyc_production = 1990
    venue = "Metropolitan Opera"
    bass_role = "Assur"
    singer = "Samuel Ramey"
    
    print(f"The bass role of {bass_role} in the {year_of_nyc_production} production of '{opera}' at the {venue} was sung by:")
    print(singer)

find_bass_singer()