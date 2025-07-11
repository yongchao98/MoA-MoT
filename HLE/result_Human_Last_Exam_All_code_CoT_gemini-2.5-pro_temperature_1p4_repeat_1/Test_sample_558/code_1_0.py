def solve_opera_puzzle():
    """
    This function solves a multi-step opera history puzzle.
    """
    # Step 1: Identify the opera based on the clues about Alboni and Caruso.
    opera = "La favorita"
    composer = "Gaetano Donizetti"
    alboni_role = "Leonora di Gusman"
    caruso_role = "Fernando"
    
    # Step 2: Determine the timeline. The Met's last staged production of this opera was in 1905.
    last_met_staged_year = 1905
    time_gap = 70
    
    # The new production was staged "more than 70 years after".
    # 1905 + 70 = 1975. So we look for a production after 1975.
    
    # Step 3: Identify the specific New York City production that fits the criteria.
    # The Opera Orchestra of New York (OONY) presented the opera in 2003.
    nyc_production_year = 2003
    nyc_production_venue = "Carnegie Hall"
    nyc_production_organization = "Opera Orchestra of New York"
    
    # Verify the timeline: 2003 - 1905 = 98 years, which is "more than 70 years".
    
    # Step 4: Identify the singer of the principal bass role in that production.
    # The main bass role in "La favorita" is Baldassarre.
    bass_role = "Baldassarre"
    bass_singer = "Vitalij Kowaljow"
    
    print(f"The opera is {opera} by {composer}.")
    print(f"The last Met-staged production with Caruso was in {last_met_staged_year}.")
    print(f"The next major NYC performance was in {nyc_production_year}, which is {nyc_production_year - last_met_staged_year} years later.")
    print(f"In the {nyc_production_year} performance by {nyc_production_organization}, the bass role of {bass_role} was sung by:")
    print(f"{bass_singer}")

solve_opera_puzzle()
<<<Vitalij Kowaljow>>>