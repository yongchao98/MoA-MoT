# Define the sets based on the premises
whales = set(['whale1', 'whale2', 'whale3'])  # All whales are artists
artists = set(['whale1', 'whale2', 'whale3', 'artist1', 'artist2'])  # Some artists are students
students = set(['artist1', 'student1'])  # Some artists are students

# Check if there is any overlap between whales and students
whales_as_students = whales.intersection(students)

# Output the result
print("Some whales are students:", len(whales_as_students) > 0)