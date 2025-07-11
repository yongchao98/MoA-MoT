import sys

# Set a higher limit for integer string conversion if necessary, although Python 3 handles this automatically.
# sys.set_int_max_str_digits(0) 

# Step 1: Assign the value for Y1 deciphered from the first clue.
# The clue points to the year of the first commercial US oil well by Edwin Drake.
Y1 = 1859

# Step 2: Assign the value for Y4 deciphered from the second clue.
# The clue is a riddle pointing to a phone number, which parses to Y4 = 3535.
Y4 = 3535

# Step 3: Perform the calculation "Y4 to the Y1".
# Python's integers have arbitrary precision, so it can handle this large exponentiation.
result = Y4 ** Y1

# Step 4: Print the final equation with all numbers, as requested.
print(f"{Y4} ^ {Y1} = {result}")

# The final answer to be extracted is the numerical result of the calculation.
final_answer_for_extraction = result