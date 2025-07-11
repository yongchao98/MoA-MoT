# This script determines and prints the hierarchy of sound importance in film dance scenes.

# The elements are ranked based on their typical role in storytelling and driving the scene's energy.
# 1. Music: Provides the core rhythm, mood, and structure.
# 2. Sound Effects: Add realism, impact, and texture (e.g., footsteps, cloth sounds).
# 3. Speech: Generally minimized or absent to avoid distracting from the music and dance.
sound_hierarchy = ["Music", "Sound Effects", "Speech"]

# Joining the list elements into a single string for clear output.
ranked_order = ", ".join(sound_hierarchy)

# Printing the final ranked order.
print(ranked_order)