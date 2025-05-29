# Total number of items for each category
num_colors = 6

# Total possible outfits
total_outfits = num_colors ** 3

# Outfits with all items the same color
same_color_outfits = num_colors

# Valid outfits
valid_outfits = total_outfits - same_color_outfits

print(valid_outfits)