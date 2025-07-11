# This script identifies the champion who can increase their hook's range
# by hooking in one direction and flashing in the opposite direction during the cast animation.

# Analysis of champions:
# - Blitzcrank, Nautilus, Pyke: Flashing repositions the origin of the hook, but does not extend its maximum travel distance in the manner described.
# - Thresh: Due to a unique interaction with his Q's wind-up animation, Flashing backward
#   before the hook is thrown will make it travel a longer distance to reach the original target point.

# Based on the analysis, only one champion fits the description.
capable_champions = ["Thresh"]

# Print the name of the champion.
# The join method is used to format the output as a comma-separated string,
# which is good practice even for a single-item list.
result = ",".join(capable_champions)
print(result)