# Extracting the numbers from the text
year_digits = [2, 0, 0, 8]
position = 5
turn = 1
act = 1

# Calculating the sum of the digits of the year
sum_of_year_digits = sum(year_digits)

# Performing the final calculation based on the derived equation
result = sum_of_year_digits - position - turn - act

# Printing the equation and the final answer
# as requested in the prompt
print(f"{year_digits[0]} + {year_digits[1]} + {year_digits[2]} + {year_digits[3]} - {position} - {turn} - {act} = {result}")