# Create a class to represent a person
class Person:
    def __init__(self, name, gender):
        self.name = name
        self.gender = gender
        self.friends = []

# Create Margaret and William
margaret = Person("Margaret", "F")
william = Person("William", "M")

# Create 6 female friends
female_friends = [Person(f"Female{i}", "F") for i in range(1, 7)]

# Connect Margaret with William
margaret.friends.append(william)
william.friends.append(margaret)

# Connect Margaret with all female friends
for female in female_friends:
    margaret.friends.append(female)
    female.friends.append(margaret)

# Connect William with all female friends
for female in female_friends:
    william.friends.append(female)
    female.friends.append(william)

# Count William's female friends
william_female_friends = sum(1 for friend in william.friends if friend.gender == "F")
print(f"William has {william_female_friends} female friends")