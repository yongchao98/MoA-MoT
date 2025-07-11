# Values of the 1992 series Belarusian ruble banknotes
bison_value = 100
elk_value = 25
lynx_value = 5
hare_value = 1

# Calculate the total value
total_value = bison_value + elk_value + lynx_value + hare_value

# Print the equation and the final sum
print(f"{bison_value} + {elk_value} + {lynx_value} + {hare_value} = {total_value}")

# The prompt also asks for just the number as the final answer.
# We will print it separately to adhere to that format if needed for parsing,
# though the above line already provides the result in context.
# In a real script, you would likely only need one of these print statements.
print(total_value)