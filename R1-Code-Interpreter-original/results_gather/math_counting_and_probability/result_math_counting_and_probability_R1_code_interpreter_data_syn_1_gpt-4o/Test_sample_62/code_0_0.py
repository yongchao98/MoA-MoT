# Initialize the string with the first digit '1'
string = '1'

# Possible two-digit numbers divisible by 19 or 31
div_19 = [19, 38, 57, 76, 95]
div_31 = [31, 62, 93]

# Convert these numbers to strings for easy concatenation
div_19_str = [str(x) for x in div_19]
div_31_str = [str(x) for x in div_31]

# Function to find the next valid digit
def find_next_digit(current_digit):
    # Check for possible continuations
    for num in div_19_str + div_31_str:
        if num.startswith(current_digit):
            return num[1]
    return None

# Construct the string
while len(string) < 2002:
    next_digit = find_next_digit(string[-1])
    if next_digit:
        string += next_digit
    else:
        break

# Output the last digit
print(string[-1])