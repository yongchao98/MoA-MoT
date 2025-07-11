# Define the properties of the ponds from the answer choices
ponds = {
    'A': {'size': 10, 'age': 1},
    'C': {'size': 30, 'age': 1},
    'D': {'size': 10, 'age': 5},
    'E': {'size': 30, 'age': 5}
}

max_threat_score = 0
best_pond_choice = ''
best_pond_details = {}

print("Calculating the medical threat score for each pond.")
print("Model: Threat Score = Surface Area (side * side) * Age (years)\n")

# Iterate through each pond to calculate its threat score
for choice, properties in ponds.items():
    side_length = properties['size']
    age = properties['age']
    
    # Calculate surface area
    surface_area = side_length * side_length
    
    # Calculate threat score
    threat_score = surface_area * age
    
    print(f"Pond {choice}:")
    print(f"  Calculation: {side_length} ft * {side_length} ft * {age} year(s) = {threat_score}")
    
    # Check if this pond has the highest score so far
    if threat_score > max_threat_score:
        max_threat_score = threat_score
        best_pond_choice = choice
        best_pond_details = properties

print("\n---")
print("Conclusion:")
print(f"The pond with the highest threat score is choice {best_pond_choice}.")
print(f"This pond is {best_pond_details['size']} feet square and {best_pond_details['age']} years old because a larger and older pond can support a larger mosquito population.")
