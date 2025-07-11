# 1. Decode the numerals into modern numbers
num1 = 4  # Represents the Chinese numeral '四'
num2 = 30 # Represents the Cuneiform numeral '𒌋𒌋𒌋' (10+10+10)
num3 = 23 # Represents the Cuneiform numeral '𒌋𒌋𒁹𒁹𒁹' (10+10+1+1+1)

# 2. Perform the calculation
result = num1 / num2 - num3

# 3. Print the equation with modern numbers and the final answer
print(f"The equation can be translated to: {num1} / {num2} - {num3}")
print(f"The result is: {result}")