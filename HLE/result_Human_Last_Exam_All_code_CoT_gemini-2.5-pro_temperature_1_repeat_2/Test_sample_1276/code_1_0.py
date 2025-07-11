# The key number derived from the flavor text is 8.
# We select the riddles where the answer length is 8.
riddle_numbers = [3, 12, 18]

# The final equation is the sum of these riddle numbers.
num1 = riddle_numbers[0]
num2 = riddle_numbers[1]
num3 = riddle_numbers[2]
result = sum(riddle_numbers)

# Print the equation
print(f"The final equation is: {num1} + {num2} + {num3} = {result}")

# The result, 33, points to the final answer phrase.
# The number 8 from the flavor text gives the first word.
final_answer = "EIGHT IS ENOUGH"
print(f"The final answer is: {final_answer}")