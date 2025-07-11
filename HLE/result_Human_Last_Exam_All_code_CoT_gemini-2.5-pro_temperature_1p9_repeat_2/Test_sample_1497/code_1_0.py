# This script determines and prints the solution based on the physical analysis
# of how different playing variations affect the parameters of a violin waveform model.

# The variations and the corresponding parameter groups they primarily affect are:
# 1. sul ponticello -> group ii (nu), related to string harmonic damping.
# 2. bridge mute -> group iv (mu, a_2, f_2), related to higher body resonances.
#    A mute adds mass, which lowers the resonance frequency f_2 ('down').
# 3. helium-filled room -> group iii (a_1, f_1), related to the main air resonance.
#    Helium increases the speed of sound, shifting the air resonance frequency f_1.
# 4. on the E string -> group i (F), the fundamental frequency of the open string.

solution = ["ii", "iv", "iii", "i", "down"]

# Print the final answer in the specified comma-separated format.
print(",".join(solution))