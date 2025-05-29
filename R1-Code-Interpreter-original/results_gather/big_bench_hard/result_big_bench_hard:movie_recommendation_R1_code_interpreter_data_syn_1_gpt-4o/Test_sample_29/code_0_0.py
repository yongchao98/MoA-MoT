# Define a dictionary with movie options and their genres/themes
movies = {
    "The Bank Job": ["Crime", "Drama", "Thriller"],
    "Robot Carnival": ["Animation", "Fantasy", "Sci-Fi"],
    "Dances with Wolves": ["Adventure", "Drama", "Western"],
    "The Family Stone": ["Comedy", "Drama", "Romance"]
}

# Define the genres/themes of the given movies
target_genres = ["Drama", "Historical", "Epic", "Crime"]

# Function to find the most similar movie
def find_similar_movie(movies, target_genres):
    similarity_scores = {}
    for movie, genres in movies.items():
        # Calculate similarity score based on common genres/themes
        similarity_scores[movie] = len(set(genres) & set(target_genres))
    # Find the movie with the highest similarity score
    similar_movie = max(similarity_scores, key=similarity_scores.get)
    return similar_movie

# Find and print the most similar movie
similar_movie = find_similar_movie(movies, target_genres)
print(similar_movie)