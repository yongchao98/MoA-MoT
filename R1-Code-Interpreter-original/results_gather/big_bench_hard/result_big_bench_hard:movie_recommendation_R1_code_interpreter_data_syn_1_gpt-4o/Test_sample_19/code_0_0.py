# Define the genres for each movie
movies = {
    "Beauty and the Beast": ["romantic", "fantasy", "musical"],
    "Independence Day": ["science fiction", "action"],
    "Forrest Gump": ["drama", "romance", "historical"],
    "Jurassic Park": ["science fiction", "adventure"]
}

# Define the genres for each option
options = {
    "The Bag Man": ["crime", "thriller"],
    "The Silence of the Lambs": ["psychological", "horror", "thriller"],
    "Extreme Measures": ["thriller", "drama"],
    "Kicking Off": ["comedy"]
}

# Find the common genres
common_genres = set(movies["Beauty and the Beast"]) | set(movies["Independence Day"]) | set(movies["Forrest Gump"]) | set(movies["Jurassic Park"])

# Calculate the match score for each option
scores = {}
for option, genres in options.items():
    match_score = len(set(genres) & common_genres)
    scores[option] = match_score

# Find the option with the highest match score
best_match = max(scores, key=scores.get)
print(best_match)