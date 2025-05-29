# Define a dictionary with movie options and their genres/themes
movies = {
    "Brave": ["Animation", "Fantasy", "Adventure"],
    "Dredd": ["Science Fiction", "Action", "Dystopian"],
    "When Worlds Collide": ["Science Fiction", "Disaster"],
    "Go": ["Crime", "Comedy", "Drama"]
}

# Define the genres/themes of the given movies
given_movies_genres = ["Historical", "Drama", "Crime", "War", "Science Fiction", "Dystopian"]

# Function to find the most similar movie
def find_similar_movie(movies, given_genres):
    similarity_scores = {}
    for movie, genres in movies.items():
        # Calculate similarity score based on common genres/themes
        similarity_scores[movie] = len(set(genres) & set(given_genres))
    # Find the movie with the highest similarity score
    similar_movie = max(similarity_scores, key=similarity_scores.get)
    return similar_movie

# Find and print the most similar movie
similar_movie = find_similar_movie(movies, given_movies_genres)
print(similar_movie)