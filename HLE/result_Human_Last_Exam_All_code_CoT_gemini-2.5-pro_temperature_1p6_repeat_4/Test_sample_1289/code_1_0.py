# Plan:
# 1. Define the 7 intervals for the ascending Maqam Zanjaran.
#    - Jins Rast (on degrees 1-4): {1, 0.75, 0.75}
#    - Connecting interval (between 4 and 5): {1}
#    - Jins Hijaz (on degrees 5-8): {0.5, 1.5, 0.5}
# 2. Define the 4 intervals for the descending path.
#    - The path is from the 8th degree down to the 4th.
#    - The scale for descent uses Jins Nahawand on the 4th degree.
#    - Jins Nahawand has intervals {1, 0.5, 1}.
#    - Let the scale degrees be: C(1), D, Eb-, F(4), G(5), Ab(6), Bb(7), C'(8).
#    - The descending path is C' -> Bb -> Ab -> G -> F.
#    - Interval 8 (C' -> Bb): 1.0
#    - Interval 9 (Bb -> Ab): 1.0
#    - Interval 10 (Ab -> G): 0.5
#    - Interval 11 (G -> F): 1.0
# 3. Combine the lists and print in the specified format.

# Ascending intervals from Maqam Zanjaran
ascending_intervals = [1, 0.75, 0.75, 1, 0.5, 1.5, 0.5]

# Descending intervals based on the modified scale
descending_intervals = [1, 1, 0.5, 1]

# Combine all 11 intervals in order
all_intervals = ascending_intervals + descending_intervals

# Format the list as a string with explicit numbers and print
# The map(str, ...) converts each number in the list to its string representation.
# ",".join(...) concatenates these strings with a comma in between.
intervals_string = ",".join(map(str, all_intervals))
print(f"{{{intervals_string}}}")
