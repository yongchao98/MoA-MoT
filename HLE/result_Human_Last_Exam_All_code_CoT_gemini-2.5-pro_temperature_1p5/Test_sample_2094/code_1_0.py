# The task is to identify which champion from the list [Blitzcrank, Nautilus, Pyke, Thresh]
# can perform a specific mechanic: hooking forward and flashing backward to extend the hook's range.

# Let's analyze each champion's ability interaction with Flash.
# Thresh (Q - Death Sentence): When Thresh casts his Q, there is a wind-up animation.
# If he Flashes during this animation, the hook originates from his new (post-Flash) location
# but travels in the direction he initially targeted.
# This means flashing backward while aiming forward creates an unexpectedly long hook.
# This matches the description perfectly.
can_perform_trick_thresh = True

# Blitzcrank (Q - Rocket Grab): When Blitzcrank Flashes during his Q animation,
# the hook's origin point moves to his new location.
# Flashing backward while aiming forward would make the hook originate from further away,
# thus DECREASING its range.
can_perform_trick_blitzcrank = False

# Nautilus (Q - Dredge Line): Works identically to Blitzcrank's hook.
# Flashing backward during the cast would decrease the hook's effective range.
can_perform_trick_nautilus = False

# Pyke (Q - Bone Skewer): The common range-extending mechanic for Pyke is to
# charge Q, Flash forward, and then release Q. The mechanic as described
# (hooking then flashing backward) does not work to extend range for him.
can_perform_trick_pyke = False

# Now, we will compile the list of champions who can perform the trick.
champions = {
    "Blitzcrank": can_perform_trick_blitzcrank,
    "Nautilus": can_perform_trick_nautilus,
    "Pyke": can_perform_trick_pyke,
    "Thresh": can_perform_trick_thresh
}

result = [champion for champion, can_do_it in champions.items() if can_do_it]

# If the list is empty, the answer is 'None'. Otherwise, join the names with a comma.
if not result:
    final_answer = "None"
else:
    final_answer = ",".join(result)

print(final_answer)