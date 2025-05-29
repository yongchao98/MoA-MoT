# Create a class to represent the friendship network
class FriendshipNetwork:
    def __init__(self):
        self.sarah_male_friends = 1  # William
        self.sarah_female_friends = 6
        
    def get_william_female_friends(self):
        # Since all Sarah's friends are friends with each other,
        # William is friends with all Sarah's female friends
        return self.sarah_female_friends

# Create the network and calculate
network = FriendshipNetwork()
william_female_friends = network.get_william_female_friends()
print(f"William has {william_female_friends} female friends")