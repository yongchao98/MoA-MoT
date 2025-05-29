# Manually constructing a possible path
def manual_word_ladder():
    # A possible path from HIYA to NEED
    path = ["HIYA", "HIDE", "HEED", "NEED"]
    return path

# Get the word ladder
ladder = manual_word_ladder()
print(",".join(ladder))