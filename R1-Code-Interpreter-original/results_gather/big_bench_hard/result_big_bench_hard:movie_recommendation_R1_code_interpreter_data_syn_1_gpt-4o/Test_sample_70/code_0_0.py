# Define the movies and their attributes
movies = {
    "American Beauty": {"genre": ["drama", "comedy"], "themes": ["suburban life", "midlife crisis"]},
    "Star Wars IV": {"genre": ["sci-fi", "adventure"], "themes": ["good vs evil", "space opera"]},
    "Star Wars V": {"genre": ["sci-fi", "adventure"], "themes": ["good vs evil", "space opera"]},
    "Forrest Gump": {"genre": ["drama", "comedy"], "themes": ["historical events", "personal growth"]},
}

options = {
    "Cargo": {"genre": ["thriller", "horror"], "themes": ["post-apocalyptic"]},
    "Neon Bull": {"genre": ["drama"], "themes": ["character study"]},
    "Jeffrey": {"genre": ["comedy", "drama"], "themes": ["love", "relationships"]},
    "The Silence of the Lambs": {"genre": ["thriller", "horror"], "themes": ["psychological"]},
}

# Function to calculate similarity score
def calculate_similarity(movie_attributes, option_attributes):
    genre_score = len(set(movie_attributes["genre"]).intersection(option_attributes["genre"]))
    theme_score = len(set(movie_attributes["themes"]).intersection(option_attributes["themes"]))
    return genre_score + theme_score

# Calculate scores for each option
scores = {option: 0 for option in options}
for movie in movies:
    for option in options:
        scores[option] += calculate_similarity(movies[movie], options[option])

# Find the option with the highest score
most_similar_movie = max(scores, key=scores.get)
print(most_similar_movie)