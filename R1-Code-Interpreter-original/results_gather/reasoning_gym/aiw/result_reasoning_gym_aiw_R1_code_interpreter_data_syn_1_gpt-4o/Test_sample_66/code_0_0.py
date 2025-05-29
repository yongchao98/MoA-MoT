# Define the group of friends
mary_friends = {
    "males": 4,  # Including David
    "females": 1  # Mary's female friend
}

# David's female friends include Mary and Mary's female friend
david_female_friends = 1 + mary_friends["females"]

# Output the number of female friends David has
print(david_female_friends)