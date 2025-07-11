# The provided image is a well-known Piet program that outputs "Hello, world!" followed by a newline.
# The following Python code constructs this output from the corresponding ASCII character codes,
# displaying the numbers involved as part of a final expression.

# The ASCII codes for "Hello, world!\n"
codes = [72, 101, 108, 108, 111, 44, 32, 119, 111, 114, 108, 100, 10]

# Build the string parts for the equation
equation_parts = [f"chr({c})" for c in codes]

# Build the final string result
result = "".join([chr(c) for c in codes])

# Print the final equation, showing all the numbers and the resulting string.
# The result string contains a newline, so we represent it as "\n" for clarity in the output.
print(f"{' + '.join(equation_parts)} = \"{result.strip()}\\n\"")

# The actual output of the Piet program is the string itself.
# We print it here for completeness.
print("\n--- Program Output ---")
print(result, end="")