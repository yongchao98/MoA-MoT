# Define a dictionary with movies and their genres/themes
movies = {
    "Pulp Fiction": ["Crime", "Drama", "Dark Humor"],
    "The Shawshank Redemption": ["Drama", "Prison", "Friendship"],
    "Dances with Wolves": ["Historical", "Drama", "Cultural Exchange"],
    "Stargate": ["Science Fiction", "Adventure"],
    "The Fugitive": ["Thriller", "Drama", "Suspense"],
    "The Boss of It All": ["Comedy"],
    "Barb Wire": ["Science Fiction", "Action"],
    "Basic Instinct 2": ["Thriller", "Psychological"]
}

# Function to find the most similar movie
def find_similar_movie(target_movies, options):
    target_genres = set()
    for movie in target_movies:
        target_genres.update(movies[movie])
    
    best_match = None
    max_common_genres = 0
    
    for option in options:
        common_genres = len(target_genres.intersection(movies[option]))
        if common_genres > max_common_genres:
            max_common_genres = common_genres
            best_match = option
    
    return best_match

# Define the target movies and options
target_movies = ["Pulp Fiction", "The Shawshank Redemption", "Dances with Wolves", "Stargate"]
options = ["The Fugitive", "The Boss of It All", "Barb Wire", "Basic Instinct 2"]

# Find and print the most similar movie
similar_movie = find_similar_movie(target_movies, options)
print(similar_movie)