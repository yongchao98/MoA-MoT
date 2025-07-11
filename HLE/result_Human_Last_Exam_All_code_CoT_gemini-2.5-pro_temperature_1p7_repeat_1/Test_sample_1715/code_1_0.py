# Step 1: Define the initial ownership state in 2016
# We use dictionaries to store the information for each person.
tommy_initial = {
    "Owner": "Tommy",
    "Lot": "A",
    "Land_Number": "1234567890",
    "Trailer_Number": "1234567890.1"
}

james_initial = {
    "Owner": "James",
    "Lot": "B",
    "Land_Number": "0987654321",
    "Trailer_Number": "0987654321.1"
}

# Step 2: Model the legal changes in 2017.
# Tommy's Estate re-records the deed to correct the lot from A to B.
tommy_final = tommy_initial.copy()
tommy_final["Owner"] = "Tommy's Estate"
tommy_final["Lot"] = "B" # The deed was corrected to Lot B.

# James was advised to correct his deed to Lot A but did not.
# However, the correction on Tommy's deed legally implies James owns the other lot (Lot A).
james_final = james_initial.copy()
james_final["Lot"] = "A" # This is the logically resulting ownership.

# Step 3: Print the final ownership status.
print("Final Ownership Status:\n")

# Print Tommy's Estate's information
print(f"Owner: {tommy_final['Owner']}")
print(f"Owns: Lot {tommy_final['Lot']}")
print(f"Land Tax Map Number: {tommy_final['Land_Number']}")
print(f"Trailer Number: {tommy_final['Trailer_Number']}\n")

# Print James's information
print(f"Owner: {james_final['Owner']}")
print(f"Owns: Lot {james_final['Lot']}")
print(f"Land Tax Map Number: {james_final['Land_Number']}")
print(f"Trailer Number: {james_final['Trailer_Number']}")

final_answer = (
    f"Tommy's Estate owns Lot B with Land Number {tommy_final['Land_Number']} and Trailer Number {tommy_final['Trailer_Number']}. "
    f"James owns Lot A with Land Number {james_final['Land_Number']} and Trailer Number {james_final['Trailer_Number']}."
)
# The output is formatted into the final answer string below.