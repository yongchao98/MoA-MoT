# Define the houses and their characteristics
houses = [
    {"name": "Alice", "animal": "bird", "smoothie": "Dragonfruit", "drink": "water"},
    {"name": "Carol", "animal": "dog", "smoothie": "Butterscotch", "drink": "milk"},
    {"name": "Bob", "animal": "cat", "smoothie": "Desert", "drink": "coffee"},
    {"name": "Arnold", "animal": "fish", "smoothie": "Darkness", "drink": "boba tea"}
]

# Check the conditions
def check_conditions(houses):
    # Condition 1: The fish enthusiast and the Butterscotch smoothie drinker are next to each other.
    condition1 = (houses[1]["smoothie"] == "Butterscotch" and houses[2]["animal"] == "fish") or \
                 (houses[2]["smoothie"] == "Butterscotch" and houses[3]["animal"] == "fish")
    
    # Condition 2: Arnold is in the fourth house.
    condition2 = houses[3]["name"] == "Arnold"
    
    # Condition 3: The cat lover is the Desert smoothie lover.
    condition3 = any(house["animal"] == "cat" and house["smoothie"] == "Desert" for house in houses)
    
    # Condition 4: The Dragonfruit smoothie lover is the one who only drinks water.
    condition4 = any(house["smoothie"] == "Dragonfruit" and house["drink"] == "water" for house in houses)
    
    # Condition 5: The Darkness smoothie drinker is the fish enthusiast.
    condition5 = any(house["smoothie"] == "Darkness" and house["animal"] == "fish" for house in houses)
    
    # Condition 6: The person who likes milk is the dog owner.
    condition6 = any(house["drink"] == "milk" and house["animal"] == "dog" for house in houses)
    
    # Condition 7: Bob is the coffee drinker.
    condition7 = any(house["name"] == "Bob" and house["drink"] == "coffee" for house in houses)
    
    # Condition 8: The cat lover is directly left of the fish enthusiast.
    condition8 = (houses[2]["animal"] == "cat" and houses[3]["animal"] == "fish")
    
    # Condition 9: Carol is in the second house.
    condition9 = houses[1]["name"] == "Carol"
    
    # Condition 10: The cat lover is Bob.
    condition10 = any(house["animal"] == "cat" and house["name"] == "Bob" for house in houses)
    
    return all([condition1, condition2, condition3, condition4, condition5, condition6, condition7, condition8, condition9, condition10])

# Check if all conditions are satisfied
print(check_conditions(houses))