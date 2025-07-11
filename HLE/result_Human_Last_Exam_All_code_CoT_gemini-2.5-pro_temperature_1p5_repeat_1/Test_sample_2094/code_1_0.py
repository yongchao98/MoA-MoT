# In League of Legends, the term "hook-flash" can mean different things.
# The user's description specifies "hooking in one direction and flashing in the opposite direction to increase the hook's range."
# This describes a specific mechanic where the ability's origin is locked before the flash.

# Let's analyze the candidates:
# Blitzcrank, Nautilus, Pyke: Their hook abilities, when combined with Flash, will originate from their location AFTER the flash.
# This means they cannot perform the specific maneuver described.

# Thresh: His hook (Death Sentence - Q) has a unique property. During the cast time, he can Flash,
# but the hook will still originate from his pre-Flash location. This allows him to cast the hook
# and then Flash backward, making the hook travel from its original point, which matches the user's description.

# Therefore, only Thresh fits the criteria.
champions_that_can_perform_the_mechanic = ["Thresh"]

# Prepare the final answer string.
answer = ",".join(champions_that_can_perform_the_mechanic)

# Print the final result.
print(answer)