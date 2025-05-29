# Define the characteristics of the given films
films = {
    "Raiders of the Lost Ark": ["action", "adventure", "historical", "strong lead"],
    "The Shawshank Redemption": ["drama", "hope", "friendship", "redemption"],
    "Inception": ["sci-fi", "thriller", "complex", "heist"],
    "Pulp Fiction": ["crime", "non-linear", "dark humor", "character development"]
}

# Define the characteristics of the options
options = {
    "Beyond the Poseidon Adventure": ["disaster", "adventure"],
    "The Chorus": ["drama", "redemption", "hope"],
    "Forrest Gump": ["drama", "adventure", "historical", "strong narrative"],
    "Scouts Guide to the Zombie Apocalypse": ["comedy", "horror"]
}

# Function to find the most similar movie
def find_similar_movie(films, options):
    # Calculate similarity based on common characteristics
    similarity_scores = {}
    for option, characteristics in options.items():
        score = 0
        for film, film_characteristics in films.items():
            score += len(set(characteristics) & set(film_characteristics))
        similarity_scores[option] = score
    # Find the option with the highest similarity score
    most_similar = max(similarity_scores, key=similarity_scores.get)
    return most_similar

# Find and print the most similar movie
most_similar_movie = find_similar_movie(films, options)
print(most_similar_movie)