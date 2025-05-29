# Original list of strings
numbers = ['-18.87', '-73', '97.04', '27.60', '25.8210', '42.6258']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [num for num in numbers if float(num) <= -22]

# Print the result
print(filtered_numbers)