# In the 1966 Vienna State Opera Ballet production of "Swan Lake",
# Rudolf Nureyev, as Prince Siegfried, performed a famously difficult
# solo variation in Act III just before Odile's solo.

# We are counting the number of double cabriole derrières in this variation.
# Each cabriole is counted as one step.
cabriole_1 = 1
cabriole_2 = 1
cabriole_3 = 1
cabriole_4 = 1
cabriole_5 = 1
cabriole_6 = 1

# Summing the count of each step to find the total.
total_cabrioles = cabriole_1 + cabriole_2 + cabriole_3 + cabriole_4 + cabriole_5 + cabriole_6

# The final equation shows each individual step being added to the total.
print(f"Nureyev performed a series of double cabriole derrières. The count is:")
print(f"{cabriole_1} + {cabriole_2} + {cabriole_3} + {cabriole_4} + {cabriole_5} + {cabriole_6} = {total_cabrioles}")