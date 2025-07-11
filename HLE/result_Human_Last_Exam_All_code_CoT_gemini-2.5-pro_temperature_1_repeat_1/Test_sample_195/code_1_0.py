# The goal is to find the least dangerous path from the adventurer '@' to the gold 'g'.
#
# 1. The adventurer '@' starts in the room on the right.
# 2. The gold 'g' is in the room on the left.
# 3. The dragon 'D' is in a room on the lower-right, representing the main danger.
#
# To find the "least dangerous" path, we must avoid the dragon.
# The adventurer must exit their current room, navigate a central corridor, and enter the gold's room.
#
# Path Analysis:
# - The adventurer's room has one exit: the door '+' below the adventurer. So the first move must be DOWN.
# - After moving down and exiting the room, the adventurer is in a long corridor. The next move must be LEFT to travel towards the gold's room.
# - The gold's room has a door on its right side. To enter it from the corridor, the adventurer must move UP to align with the door.
# - Finally, after moving up, the adventurer can move LEFT through the door to reach the gold.
#
# This path (Down, Left, Up, Left) keeps the adventurer in the upper part of the map,
# far away from the dragon located in the lower section. This makes it the least dangerous route.

path = "DLUL"
print(path)
