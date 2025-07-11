import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# This script deduces the 16-character string that identifies the parameters for each plot.
# The logic is based on analyzing the visual features of each plot and mapping them
# to the physical effects of the parameters b, c, and d in the nonlinear wave equation.

# Final mapping based on detailed analysis:
# P1: Z. (1,0,1). b=1,d=1 => reddish bias.
# P2: C. (-1,1,-1). c=1=>chaotic, b=-1,d=-1 => blueish.
# P3: z. (0,1,1). c=1=>chaotic, d=1=>reddish. b=0 unique.
# P4: c. (0,-1,0). c=-1=>wavy, b=0,d=0 => symmetric.
# P5: z. (0,-1,-1). c=-1=>wavy, d=-1=>blue bg. b=0 unique.
# P6: 0. (-1,-1,0). c=-1=>wavy, b=-1=>blueish. d=0 unique.
# P7: 0. (1,1,0). c=1=>chaotic, b=1=>reddish. d=0 unique.
# P8: c. (1,-1,1). c=-1=>wavy, b=1,d=1=>reddish.
# P9: B. (1,0,0). b=1=>reddish.
# P10: b. (-1,0,0). b=-1=>blueish.
# P11: C. (0,1,0). c=1=>chaotic, symmetric.
# P12: d. (0,0,-1). d=-1=>blue bg.
# P13: Z. (-1,0,-1). b=-1,d=-1=>very blue bg. c=0 unique.
# P14: B. (1,-1,-1). c=-1=>wavy, b=1,-d=-1=>complex pattern.
# P15: D. (0,0,1). d=1=>positive yellow bg.
# P16: D. (-1,-1,1). c=-1=>wavy, b=-1 vs d=1=>wavy blueish pattern.

solution_map = {
    1: 'Z',
    2: 'C',
    3: 'z',
    4: 'c',
    5: 'z',
    6: '0',
    7: '0',
    8: 'c',
    9: 'B',
    10: 'b',
    11: 'C',
    12: 'd',
    13: 'Z',
    14: 'B',
    15: 'D',
    16: 'D'
}

# Constructing the final string in order from plot 1 to 16.
final_string = ""
for i in range(1, 17):
    final_string += solution_map[i]

print(f"{final_string}")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Final print to the user
print(output)