# Define a function to check the genre and cultural impact of the movies
def check_movie_similarity():
    # Dictionary of movies with their genres and cultural impact
    movies = {
        "American Beauty": {"genre": ["Drama"], "cultural_impact": "High"},
        "Austin Powers: International Man of Mystery": {"genre": ["Comedy"], "cultural_impact": "High"},
        "Star Wars Episode V - The Empire Strikes Back": {"genre": ["Science Fiction"], "cultural_impact": "High"},
        "Mission Impossible": {"genre": ["Action", "Spy"], "cultural_impact": "High"},
        "Dragon Ball The Curse Of The Blood Rubies": {"genre": ["Animated", "Fantasy"], "cultural_impact": "Medium"},
        "Moonlight Mile": {"genre": ["Drama"], "cultural_impact": "Low"},
        "Hollywood Homicide": {"genre": ["Comedy", "Action"], "cultural_impact": "Medium"},
        "Star Wars Episode IV - A New Hope": {"genre": ["Science Fiction"], "cultural_impact": "High"}
    }
    
    # List of original movies
    original_movies = ["American Beauty", "Austin Powers: International Man of Mystery", 
                       "Star Wars Episode V - The Empire Strikes Back", "Mission Impossible"]
    
    # Calculate the similarity score for each option
    similarity_scores = {}
    for option in ["Dragon Ball The Curse Of The Blood Rubies", "Moonlight Mile", 
                   "Hollywood Homicide", "Star Wars Episode IV - A New Hope"]:
        score = 0
        for movie in original_movies:
            # Check genre overlap
            if any(genre in movies[option]["genre"] for genre in movies[movie]["genre"]):
                score += 1
            # Check cultural impact
            if movies[option]["cultural_impact"] == movies[movie]["cultural_impact"]:
                score += 1
        similarity_scores[option] = score
    
    # Find the option with the highest similarity score
    most_similar_movie = max(similarity_scores, key=similarity_scores.get)
    return most_similar_movie

# Print the most similar movie
print(check_movie_similarity())