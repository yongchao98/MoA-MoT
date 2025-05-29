# Define a dictionary with movie options and their genres
movie_genres = {
    "Care Bears Movie II A New Generation": ["Animation", "Fantasy", "Family"],
    "Noises Off": ["Comedy"],
    "The Fugitive": ["Action", "Thriller", "Crime"],
    "Lights Out": ["Horror", "Thriller"]
}

# Define the common genres we are looking for
common_genres = ["Action", "Thriller", "Crime"]

# Function to find the most similar movie based on genres
def find_similar_movie(movie_genres, common_genres):
    for movie, genres in movie_genres.items():
        if any(genre in genres for genre in common_genres):
            return movie

# Find and print the most similar movie
similar_movie = find_similar_movie(movie_genres, common_genres)
print(similar_movie)