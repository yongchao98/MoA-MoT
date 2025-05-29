# Create a class to represent the friendship network
class FriendshipNetwork:
    def __init__(self):
        self.network = {}
    
    def add_person(self, name, gender):
        self.network[name] = {'gender': gender, 'friends': set()}
    
    def add_friendship(self, person1, person2):
        self.network[person1]['friends'].add(person2)
        self.network[person2]['friends'].add(person1)
    
    def get_female_friends(self, person):
        return sum(1 for friend in self.network[person]['friends'] 
                  if self.network[friend]['gender'] == 'F')

# Create the network
network = FriendshipNetwork()

# Add people
network.add_person('Susan', 'F')
network.add_person('John', 'M')
network.add_person('Friend1', 'F')
network.add_person('Friend2', 'F')
network.add_person('Friend3', 'F')

# Add all friendships (everyone is friends with everyone)
people = ['Susan', 'John', 'Friend1', 'Friend2', 'Friend3']
for i in range(len(people)):
    for j in range(i+1, len(people)):
        network.add_friendship(people[i], people[j])

# Count John's female friends
john_female_friends = network.get_female_friends('John')
print(f"John has {john_female_friends} female friends")