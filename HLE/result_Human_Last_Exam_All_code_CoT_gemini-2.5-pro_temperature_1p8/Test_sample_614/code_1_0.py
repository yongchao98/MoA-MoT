# Let's solve the riddle step by step to find each number's position.
# "Number 2 is the best of them;" -> The first number is 2.
first_number = 2

# "Number 4 likes two the most and it always goes behind." -> The second number is 4.
second_number = 4

# "...the fifth, who never had a thing and lastly follows." -> The last number is 5.
fifth_number = 5

# "Number three is the strongest and it protects the last;" -> The fourth number is 3.
fourth_number = 3

# "Number 1 is the first to have something but is the weakest;"
# By elimination, the remaining number 1 fills the only empty spot, the third position.
third_number = 1

# Print the final sequence of numbers.
print(f"The correct order is:")
print(f"{first_number} {second_number} {third_number} {fourth_number} {fifth_number}")