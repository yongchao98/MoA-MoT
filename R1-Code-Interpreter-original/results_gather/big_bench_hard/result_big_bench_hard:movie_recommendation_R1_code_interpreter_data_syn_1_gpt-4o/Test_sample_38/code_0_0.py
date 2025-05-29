# Define the genres and themes for each movie
movies = {
    "Pulp Fiction": ["crime", "drama", "dark comedy"],
    "The Shawshank Redemption": ["drama", "hope", "friendship", "redemption"],
    "Aladdin": ["animated", "musical", "fantasy", "adventure", "romance"],
    "The Lion King": ["animated", "musical", "drama", "family", "coming of age"],
}

options = {
    "Terminator 2 Judgment Day": ["science fiction", "action", "technology", "survival"],
    "The Next Three Days": ["thriller", "crime", "drama"],
    "Detachment": ["drama", "education", "personal struggle"],
    "Ruby Sparks": ["romantic comedy", "drama", "fantasy"],
}

# Calculate similarity based on common genres/themes
def calculate_similarity(movie_themes, option_themes):
    return len(set(movie_themes) & set(option_themes))

# Calculate total similarity for each option
similarity_scores = {}
for option, themes in options.items():
    total_similarity = sum(calculate_similarity(themes, movie_themes) for movie_themes in movies.values())
    similarity_scores[option] = total_similarity

# Find the option with the highest similarity score
most_similar_movie = max(similarity_scores, key=similarity_scores.get)
print(most_similar_movie)