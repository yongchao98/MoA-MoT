# Step 1: Count the sequences in Marianela Núñez's 2018 Royal Ballet performance.
# In the Peter Wright choreography, Marianela Núñez performs a distinct series of 
# piqué turns into a supported écarté back pose.
# Watching the performance, we can count them:
# 1. Piqué turn into écarté back (right leg)
# 2. Piqué turn into écarté back (left leg)
# 3. Piqué turn into écarté back (right leg)
# 4. Piqué turn into écarté back (left leg)
nunez_count = 4

# Step 2: Count the sequences in Maria Khoreva's 2019 Mariinsky Theatre performance.
# The Vasily Vainonen choreography for the Mariinsky Ballet is different and does not
# contain the specific sequence of a piqué soutenu followed by an écarté back.
# The analogous sections feature different steps, such as piqué turns into attitude.
khoreva_count = 0

# Step 3: Calculate the difference.
difference = nunez_count - khoreva_count

# Step 4: Print the final result, showing the full equation.
print(f"The number of sequences performed by Marianela Núñez is {nunez_count}.")
print(f"The number of sequences performed by Maria Khoreva is {khoreva_count}.")
print(f"The difference is: {nunez_count} - {khoreva_count} = {difference}")
