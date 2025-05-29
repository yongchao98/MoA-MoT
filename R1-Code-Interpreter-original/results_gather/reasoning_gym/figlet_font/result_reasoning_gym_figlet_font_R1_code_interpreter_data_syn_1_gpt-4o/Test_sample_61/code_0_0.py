# Define the list of numbers
numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# Initialize the sum of even numbers
sum_of_evens = 0

# Iterate through the list and sum the even numbers
for number in numbers:
    if number % 2 == 0:
        sum_of_evens += number

# Print the result
print(sum_of_evens)