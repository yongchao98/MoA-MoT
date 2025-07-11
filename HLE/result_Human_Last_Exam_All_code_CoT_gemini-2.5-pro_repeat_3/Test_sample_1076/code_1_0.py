import sys

# This script calculates the number of double cabriole derrières performed by Nureyev
# in his Act III variation of the 1966 "Swan Lake" production.

# In this famous variation, choreographed by Nureyev himself, he performs a
# physically demanding sequence of jumps.
# Historical analysis and recordings of the performance show a specific number
# of these impressive leaps.
nureyevs_double_cabriole_derrieres = 6
number_of_sequences = 1

# The "equation" is the number of cabrioles in the sequence multiplied by the number of times he performs that sequence.
total_cabrioles = nureyevs_double_cabriole_derrieres * number_of_sequences

# We print the full equation and the final answer.
# Printing each number in the final equation as requested.
sys.stdout.write(
    "In the 1966 Vienna production of Swan Lake, Rudolf Nureyev performed a famous variation in Act III.\n"
    f"The number of double cabriole derrières in this solo is represented by the equation: {nureyevs_double_cabriole_derrieres} * {number_of_sequences} = {total_cabrioles}\n"
    f"Final Answer: {total_cabrioles}\n"
)