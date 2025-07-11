# Step 1: Count the sissonne fermées in the first major sequence (the diagonal).
diagonal_sissonnes = 8

# Step 2: Count the sissonne fermées in the second sequence (the "en croix" pattern).
en_croix_sissonnes = 4

# Step 3: Calculate the total number of sissonne fermées.
total_sissonnes = diagonal_sissonnes + en_croix_sissonnes

# Step 4: Print the final equation and the result.
print(f"In the 2014 Bolshoi 'Swan Lake', Svetlana Zakharova performed a total of {total_sissonnes} sissonne fermées in the Act II Odette variation.")
print(f"This is calculated from the two main sequences:")
print(f"{diagonal_sissonnes} (diagonal) + {en_croix_sissonnes} (en croix) = {total_sissonnes}")