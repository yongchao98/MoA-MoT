# Initial counts of items
yellow_items = {'plates': 1, 'textbooks': 2}
turquoise_items = {'plates': 3, 'textbooks': 3}
burgundy_items = {'plates': 1, 'textbooks': 3}

# Remove all plates
yellow_items['plates'] = 0
turquoise_items['plates'] = 0
burgundy_items['plates'] = 0

# Count remaining yellow items
remaining_yellow_items = yellow_items['textbooks']

print(remaining_yellow_items)