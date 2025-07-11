import sys
# A function that serves as a knowledge base for this specific opera trivia.
def find_bass_singer():
    """
    This function solves a multi-step opera trivia question.
    It identifies the opera, the specific production, and the singer for the bass role.
    """
    
    # Step 1: Identify the opera using clues about Alboni and Caruso.
    opera_title = "Lucrezia Borgia by Gaetano Donizetti"
    alboni_role = "Maffio Orsini"
    caruso_role = "Gennaro"
    caruso_last_performance_year = 1905
    
    # Step 2: Identify the specific New York City production.
    # The clue states "more than 70 years after" Caruso's performance (1905).
    # 1905 + 70 = 1975.
    # The Metropolitan Opera's 2011 production is the first since the Caruso era.
    nyc_production_venue = "Metropolitan Opera"
    nyc_production_year = 2011
    
    # Step 3: Identify the principal bass role and the singer in that production.
    bass_role = "Don Alfonso, Duke of Ferrara"
    bass_singer = "Vitalij Kowaljow"
    
    print(f"Opera: {opera_title}")
    print(f"Caruso's last performance at the Met in this opera was in {caruso_last_performance_year}.")
    print(f"The New York City production referenced is the {nyc_production_venue} staging in {nyc_production_year}.")
    print(f"The bass role in this opera is {bass_role}.")
    print(f"The singer who performed the bass role in this production was: {bass_singer}")

find_bass_singer()