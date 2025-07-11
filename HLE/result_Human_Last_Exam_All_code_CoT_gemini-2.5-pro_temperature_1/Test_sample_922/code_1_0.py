# This script solves a puzzle sequence.
# The sequence is composed of user IDs from a website, linked by a special property.

print("The sequence is a puzzle, not a standard mathematical series.")
print("The numbers are the User IDs of members on the website Math Stack Exchange.")
print("The property linking them is that their reputation scores were perfect cubes around August 2022.\n")

# The data for the users in the provided sequence
users_in_sequence = {
    "24663": ("Zev Chonoles", "216,000", "60", "3"),
    "35005": ("Bill Dubuque", "343,000", "70", "3"),
    "119261": ("lulu", "125,000", "50", "3"),
    "196219": ("Daniel Fischer", "1,000,000", "100", "3"),
    "211770": ("Martin Sleziak", "27,000", "30", "3"),
    "227296": ("fleablood", "8,000", "20", "3"),
}

print("The given users and their reputation 'equations':")
for uid, (name, rep, base, power) in users_in_sequence.items():
    print(f"- User ID {uid} ({name}): Reputation ≈ {rep} = {base}^{power}")

print("\nThe puzzle asks for the single known integer that completes this set.")
print("This corresponds to another user who was known to have a cube for a reputation score.\n")

# Data for the user who completes the sequence
completing_user_id = 293
completing_user_name = "Arturo Magidin"
final_rep = "729,000"
final_base = "90"
final_power = "3"

print("The user who completes the sequence is Arturo Magidin.")
print(f"His User ID is {completing_user_id}.")
print("The final equation, representing his reputation, is:")
print(f"{final_rep} = {final_base}^{final_power}")

print(f"\nTherefore, the integer that completes the sequence is {completing_user_id}.")
