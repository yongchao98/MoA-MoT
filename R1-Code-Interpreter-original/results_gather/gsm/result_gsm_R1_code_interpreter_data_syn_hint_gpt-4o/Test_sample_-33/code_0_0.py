# Given number of red squares
red_squares = 3347396

# Calculate the number of squares for each color
blue_squares = red_squares + 4
green_squares = blue_squares + 6
white_squares = green_squares - 12

# Total number of squares
total_squares = red_squares + blue_squares + green_squares + white_squares

# Output the total number of square feet
print(total_squares)