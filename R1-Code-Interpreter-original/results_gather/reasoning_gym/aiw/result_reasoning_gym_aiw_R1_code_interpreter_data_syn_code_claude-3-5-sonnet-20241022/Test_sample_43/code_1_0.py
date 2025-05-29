# Create a class to represent the friendship network
class FriendshipNetwork:
    def __init__(self):
        # Initialize Susan and Joseph
        self.friends = {
            'Susan': {'gender': 'F', 'friends': set()},
            'Joseph': {'gender': 'M', 'friends': set()}
        }
        # Add 6 female friends (excluding Susan)
        for i in range(6):
            name = f'Female{i+1}'
            self.friends[name] = {'gender': 'F', 'friends': set()}
        
        # Connect all friends
        for person1 in self.friends:
            for person2 in self.friends:
                if person1 != person2:
                    self.friends[person1]['friends'].add(person2)
                    self.friends[person2]['friends'].add(person1)
    
    def count_female_friends(self, person):
        # Count female friends excluding Susan
        return sum(1 for friend in self.friends[person]['friends'] 
                  if self.friends[friend]['gender'] == 'F' and friend != 'Susan')

# Create network and count Joseph's female friends
network = FriendshipNetwork()
joseph_female_friends = network.count_female_friends('Joseph')
print(f"Joseph has {joseph_female_friends} female friends")