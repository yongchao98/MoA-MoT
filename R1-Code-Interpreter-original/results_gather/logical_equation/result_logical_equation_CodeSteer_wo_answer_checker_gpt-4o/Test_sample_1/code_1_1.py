# Possible values for the letters
values = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Constraints checking function
def satisfies_constraints(assignment):
    # Extract values from the assignment
    A = assignment.get('A')
    B = assignment.get('B')
    C = assignment.get('C')
    D = assignment.get('D')
    E = assignment.get('E')
    F = assignment.get('F')
    H = assignment.get('H')
    
    # Check each constraint
    if C is not None and B is not None and C != 3 * B:
        return False
    if A is not None and C is not None and A != 4 * C:
        return False
    if F is not None and A is not None and F <= A:
        return False
    if D is not None and H is not None and D != 4 * H:
        return False
    if E is not None and B is not None and E <= B:
        return False
    if A is not None and H is not None and A - H != 32:
        return False
    if F is not None and H is not None and F + H != 84:
        return False
    
    return True

# Backtracking function
def backtrack(assignment, letters, index):
    if index == len(letters):
        return assignment
    
    letter = letters[index]
    
    for value in values:
        if value not in assignment.values():
            assignment[letter] = value
            if satisfies_constraints(assignment):
                result = backtrack(assignment, letters, index + 1)
                if result is not None:
                    return result
            del assignment[letter]
    
    return None

# List of letters in alphabetical order
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

# Start backtracking
assignment = backtrack({}, letters, 0)

# Prepare the final list in alphabetical order
final_list = [assignment[letter] for letter in letters]

# Print the result
print(f"<<<{final_list}>>>")