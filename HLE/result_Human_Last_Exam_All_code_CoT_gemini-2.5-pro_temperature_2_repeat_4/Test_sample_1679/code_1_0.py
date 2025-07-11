def find_movie_with_obelisk():
    """
    This script explains the reasoning to find the first Best Picture winner
    depicting a Luxor Obelisk and prints the answer.
    """
    
    movie_title = "An American in Paris"
    movie_year = 1951
    obelisk_location = "Place de la Concorde, Paris"

    print("Step 1: Identifying the Luxor Obelisks and the scope of the search.")
    print("The two Luxor Obelisks are located in Luxor, Egypt, and in the Place de la Concorde in Paris, France.")
    print("The search for the answer involves checking the list of Academy Award winners for Best Picture for scenes depicting either of these locations.\n")
    
    print("Step 2: Analyzing early Best Picture winners set in Paris.")
    print("Several early winners are set wholly or partially in Paris, such as 'The Life of Emile Zola' (1937) and 'Casablanca' (1943).")
    print("However, a review of these films shows they do not contain a clear, prominent depiction of the Luxor Obelisk in the Place de la Concorde.\n")
    
    print("Step 3: Identifying the first winner with a clear depiction.")
    print(f"The first Best Picture winner to feature a definitive depiction is '{movie_title}', which won for the year {movie_year}.")
    print(f"In the film's celebrated 17-minute closing ballet sequence, a stylized version of the {obelisk_location} is featured.")
    print("The set, designed in the style of French painter Raoul Dufy, clearly includes the famous obelisk as a central element.\n")

    print("Conclusion:")
    print(f"The first winner of the Academy Award for Best Picture to depict a Luxor Obelisk on-screen was '{movie_title}'.")

find_movie_with_obelisk()
