# The Beyond-Use Date (BUD) for sterile preparations is determined by USP <797> guidelines.
# These guidelines have specific rules for different container types and environments.

# Define the parameters based on the user's question.
container_type = "single dose ampule"
environment = "sterile (ISO Class 5)"

# According to USP <797>, an ampule is a sealed container that cannot be resealed after opening.
# This makes it an open system, and its contents must be used immediately.
# Any unused medication must be discarded.
# Therefore, the allowable time for use after opening is effectively zero.
bud_in_hours = 0

print(f"Calculating the BUD for a '{container_type}' opened in a '{environment}':")
print("Rule: An opened ampule is an open system and must be used immediately.")
print("Any unused portion must be discarded.")
print(f"Resulting BUD from time of opening = {bud_in_hours} hours (Immediate Use)")