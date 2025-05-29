# Define the movies and their characteristics
movies = {
    "Dances with Wolves": ["historical drama", "cultural exchange", "transformation"],
    "Stargate": ["science fiction", "adventure", "exploration"],
    "Pulp Fiction": ["crime", "non-linear narrative", "dark humor"],
    "The Shawshank Redemption": ["drama", "hope", "redemption"],
    "The Towering Inferno": ["disaster", "survival", "heroism"],
    "The Collector": ["psychological thriller", "kidnapping", "suspense"],
    "Diabolique": ["psychological thriller", "murder plot", "suspense"],
    "The Fugitive": ["action thriller", "wrongful accusation", "redemption"]
}

# Function to find the most similar movie
def find_similar_movie(target_movies, options):
    target_genres = set()
    for movie in target_movies:
        target_genres.update(movies[movie])
    
    best_match = None
    max_overlap = 0
    for option in options:
        overlap = len(target_genres.intersection(movies[option]))
        if overlap > max_overlap:
            max_overlap = overlap
            best_match = option
    
    return best_match

# Define the target movies and options
target_movies = ["Dances with Wolves", "Stargate", "Pulp Fiction", "The Shawshank Redemption"]
options = ["The Towering Inferno", "The Collector", "Diabolique", "The Fugitive"]

# Find and print the most similar movie
similar_movie = find_similar_movie(target_movies, options)
print(similar_movie)