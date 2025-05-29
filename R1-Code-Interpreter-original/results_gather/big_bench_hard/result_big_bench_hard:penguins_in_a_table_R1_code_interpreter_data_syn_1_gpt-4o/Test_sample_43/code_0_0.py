# List of ages for penguins and giraffes
penguin_ages = [7, 5, 9, 8, 12]
giraffe_ages = [5, 10, 2, 9]

# Count animals older than 5 years
count = sum(age > 5 for age in penguin_ages + giraffe_ages)

print(count)