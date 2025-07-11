# The piece is Bach's Prelude in C Major, BWV 846.
bwv_number = 846

# Convert the number to a string to iterate through its digits.
s_bwv_number = str(bwv_number)

# Initialize a list to hold the digits and the sum.
digits = []
total_sum = 0

# Calculate the sum of the digits.
for digit_char in s_bwv_number:
    digit_int = int(digit_char)
    digits.append(digit_int)
    total_sum += digit_int

# Print the equation and the result.
equation = " + ".join(map(str, digits))
print(f"The piece is BWV {bwv_number}.")
print(f"Summing the digits of the BWV number: {equation} = {total_sum}")

# The result, 18, corresponds to one of the answer choices.
# This opus number, Op. 18, famously belongs to Rachmaninoff's Piano Concerto No. 2.