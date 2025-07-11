import sys
import io

# Backup stdout
stdout_backup = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = io.StringIO()

# The plan is to define the interval sequences for the ascending and descending parts
# of the musical phrase, combine them, and then format the result as requested.

# 1. Define the ascending intervals based on Maqam Zanjaran (Jins Hijaz + Jins Nikriz)
# Intervals 1-3 from Jins Hijaz on the tonic: 0.5, 1.5, 0.5
# Intervals 4-6 from Jins Nikriz on the 4th: 0.5, 1.5, 0.5
# Interval 7 from 7th to octave: 1.0
ascending_intervals = [0.5, 1.5, 0.5, 0.5, 1.5, 0.5, 1.0]

# 2. Define the descending intervals based on the modified scale (with Jins Nahawand on the 4th)
# The descent is from note 8 down to 4.
# The upper part of the scale is defined by Jins Nahawand on 4th:
# Interval 4->5 = 1.0, 5->6 = 0.5, 6->7 = 1.0. The 7->8 interval is 1.0.
# The descending intervals are:
# Note 8 -> 7: 1.0
# Note 7 -> 6: 1.0
# Note 6 -> 5: 0.5
# Note 5 -> 4: 1.0
descending_intervals = [1.0, 1.0, 0.5, 1.0]

# 3. Combine the two lists to get the full sequence of 11 intervals.
all_intervals = ascending_intervals + descending_intervals

# 4. Format and print the result as a single string.
# The format must be {n1,n2,n3,...} without any spaces.
output_string = "{" + ",".join(map(str, all_intervals)) + "}"
# We need to print this to the captured stdout
# Redirect back to original stdout
sys.stdout = stdout_backup
# Print the captured output
print(output_string)

# The final answer in the requested format
final_answer = output_string
print(f'<<<{final_answer}>>>')