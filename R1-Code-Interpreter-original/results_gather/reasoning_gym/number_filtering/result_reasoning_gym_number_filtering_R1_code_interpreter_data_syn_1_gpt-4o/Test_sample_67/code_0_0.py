# Original list of strings
original_list = ['-23.8495', '-64.2', '-73', '-49.770', '6.1625', '-82.0327', '-2.145', '44.4451', '-56.66']

# Convert strings to floats and filter out numbers greater than 22.0
filtered_list = [str(num) for num in map(float, original_list) if num <= 22.0]

# Print the filtered list
print(filtered_list)