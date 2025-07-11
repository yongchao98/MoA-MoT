# Set the dimension of the system
d = 3

# Calculate the average number of contacts for a particle in a jammed state
# In this state, the system is isostatic.
avg_contacts_jammed = 2 * d

# Calculate the number of contacts for the unstable particles
# These particles are not strictly jammed and can move.
contacts_unstable = d + 1

print(f"In a hard-sphere system at jamming in d={d} dimensions:")
print(f"The average number of contacts for a stable, jammed particle is 2 * {d} = {avg_contacts_jammed}.")
print(f"The number of contacts for the unstable particles described is {d} + 1 = {contacts_unstable}.")
