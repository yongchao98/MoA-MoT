import collections

# Assignments based on physical reasoning and visual inspection
assignments = {
    1: 'B', 2: 'b', 3: 'D', 4: 'b',
    5: 'd', 6: 'c', 7: 'D', 8: 'c',
    9: 'B', 10: 'C', 11: 'Z', 12: 'd',
    13: 'b', 14: 'Z', 15: 'C', 16: '0'
}

# Sort by plot number to ensure correct order
sorted_assignments = collections.OrderedDict(sorted(assignments.items()))

# Construct the final answer string
answer_string = "".join(sorted_assignments.values())

print(answer_string)