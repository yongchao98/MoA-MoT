# Define a function to check the genres and themes of each movie
def check_movie_similarity():
    # Define the movies and their genres/themes
    movies = {
        "Dragon Ball Z Dead Zone": ["Fantasy", "Action", "Animation"],
        "Lightning in a Bottle": ["Documentary", "Music"],
        "The Omega Man": ["Science Fiction", "Action", "Post-Apocalyptic"],
        "Goodfellas": ["Crime", "Drama", "Biography"]
    }
    
    # Define the common themes/genres of the given movies
    common_themes = ["Fantasy", "Unique Storytelling", "Cult Classic", "Surrealism"]
    
    # Check which movie has the most similar themes/genres
    similarity_scores = {}
    for movie, genres in movies.items():
        similarity_scores[movie] = len(set(genres) & set(common_themes))
    
    # Find the movie with the highest similarity score
    most_similar_movie = max(similarity_scores, key=similarity_scores.get)
    return most_similar_movie

# Print the most similar movie
print(check_movie_similarity())