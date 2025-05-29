# Define the friends
susan_friends = {
    "males": ["Michael", "MaleFriend2", "MaleFriend3"],
    "females": ["FemaleFriend1"]
}

# Michael's friends
michael_friends = susan_friends["males"] + susan_friends["females"] + ["Susan"]

# Count female friends of Michael
female_friends_of_michael = [friend for friend in michael_friends if friend in susan_friends["females"] or friend == "Susan"]

# Output the number of female friends Michael has
print(len(female_friends_of_michael))