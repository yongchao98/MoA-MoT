# Define the conditions based on the premises
def is_uncle_of_owen(person):
    return person in ["uncle_of_owen"]

def is_nephew_of_gene(person):
    return person in ["schoolmate_of_steve"]

def is_nephew_of_ignacio_and_brother_of_david(person):
    return person in ["nephew_of_ignacio_and_brother_of_david"]

# Check the conclusion
def check_conclusion(person):
    if is_nephew_of_ignacio_and_brother_of_david(person):
        # They must be either a schoolmate of Steve or an uncle of Owen
        if is_nephew_of_gene(person):
            return False  # Contradicts the conclusion
        if is_uncle_of_owen(person):
            return True  # Satisfies the conclusion
    return False

# Test the logic
person = "nephew_of_ignacio_and_brother_of_david"
result = check_conclusion(person)
print("Conclusion holds:", result)