# 1. Define the numbers based on historical and numeral system interpretation.
gwalior_symbol = " ग्वालियर"
numeral_1_symbol = "𒌋𒌋𒌋"
numeral_2_symbol = "𒌋𒌋𒁹𒁹𒁹"

# In the Gwalior numeral system and its famous inscription, 'Gwalior' is linked to the first appearance of '0' as a symbol.
gwalior_value = 0
# The other symbols are Babylonian cuneiform. '𒌋' is 10 and '𒁹' is 1.
numeral_1_value = 10 + 10 + 10
numeral_2_value = 10 + 10 + 1 + 1 + 1

# 2. Explain the conversion.
print(f"Decoding the equation: {gwalior_symbol} / {numeral_1_symbol} - {numeral_2_symbol}")
print("------------------------------------------")
print(f"'{gwalior_symbol}' represents {gwalior_value}, referencing the Gwalior inscription's early use of zero.")
print(f"'{numeral_1_symbol}' is Cuneiform for 10 + 10 + 10 = {numeral_1_value}.")
print(f"'{numeral_2_symbol}' is Cuneiform for 10 + 10 + 1 + 1 + 1 = {numeral_2_value}.")
print("------------------------------------------\n")


# 3. Formulate and solve the equation using the decoded numbers.
print("The translated equation is:")
result = gwalior_value / numeral_1_value - numeral_2_value

# 4. Print the final equation with all numbers and the result.
print(f"{gwalior_value} / {numeral_1_value} - {numeral_2_value} = {result}")

# The final answer in the requested format
print(f"<<<{-23.0}>>>")