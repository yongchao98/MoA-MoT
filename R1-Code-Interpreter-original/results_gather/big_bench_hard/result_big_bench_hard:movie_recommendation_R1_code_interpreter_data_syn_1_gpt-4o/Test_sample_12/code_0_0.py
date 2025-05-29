# Define the genres for each movie
movies = {
    "Forrest Gump": ["Drama", "Romance", "Comedy"],
    "The Silence of the Lambs": ["Thriller", "Horror"],
    "Mission Impossible": ["Action", "Thriller"],
    "Jurassic Park": ["Science Fiction", "Adventure"],
    "Joe Somebody": ["Comedy", "Drama"],
    "Dogfight": ["Romantic", "Drama"],
    "Independence Day": ["Science Fiction", "Action"],
    "Twin Peaks Fire Walk with Me": ["Horror", "Thriller"]
}

# Define the target genres based on the given movies
target_genres = set(movies["Forrest Gump"] + movies["The Silence of the Lambs"] + 
                    movies["Mission Impossible"] + movies["Jurassic Park"])

# Function to calculate similarity based on genre overlap
def calculate_similarity(movie_genres, target_genres):
    return len(set(movie_genres) & target_genres)

# Calculate similarity for each option
similarities = {movie: calculate_similarity(genres, target_genres) for movie, genres in movies.items() if movie in ["Joe Somebody", "Dogfight", "Independence Day", "Twin Peaks Fire Walk with Me"]}

# Find the movie with the highest similarity
most_similar_movie = max(similarities, key=similarities.get)
print(most_similar_movie)