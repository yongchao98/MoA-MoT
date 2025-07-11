# This script prints the solution to the violin physics problem.
# The solution is a string that identifies the parameter group changed in each of the four scenarios,
# followed by the direction of change for the parameter in the second scenario.

# The reasoning for the answer is as follows:
# 1. Sul ponticello (playing near the bridge) changes the harmonic content of the string's vibration, corresponding to group ii (ν).
# 2. A bridge mute adds mass to the bridge, lowering the body's resonant frequencies, corresponding to group iii (a1, f1). The frequency f1 goes down.
# 3. Helium in the room changes the medium in the violin body, increasing the air resonance frequency, corresponding to group iv (a2, f2).
# 4. Playing on the E string instead of the A string changes the open string's fundamental frequency, corresponding to group i (F).

# The final answer is formatted as a comma-separated string: group1,group2,group3,group4,direction
solution_string = "ii,iii,iv,i,down"

print(solution_string)