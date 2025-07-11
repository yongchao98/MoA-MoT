# Initial state of ownership in 2016
tommy_lot = "A"
tommy_number = "1234567890"
tommy_trailer_number = "1234567890.1"

james_lot = "B"
james_number = "0987654321"
james_trailer_number = "0987654321.1"

# In 2017, Tommy's Estate corrected the deed.
# Tommy's deed was re-filed and re-recorded from Lot A to Lot B.
# This is a legally binding change.
tommy_lot = "B"

# James was advised to correct his deed to Lot A, but never did.
# However, since Tommy's Estate now legally owns Lot B via the corrected deed,
# Lot A is the remaining parcel. By legal succession and the intent of the
# correction, Lot A belongs to James, even if his paperwork is outdated.
james_lot = "A"

# The final ownership status based on the legally recorded changes.
print("Final Ownership Status:")

# Output Tommy's Estate's ownership details
print(f"Tommy's Estate owns Lot {tommy_lot}")
print(f"  Land Tax Map Number: {tommy_number}")
print(f"  Trailer Tax Map Number: {tommy_trailer_number}")

print("-" * 20)

# Output James' ownership details
print(f"James owns Lot {james_lot}")
print(f"  Land Tax Map Number: {james_number}")
print(f"  Trailer Tax Map Number: {james_trailer_number}")

# Prepare the final answer in the requested format
final_answer = (
    f"Tommy's Estate owns Lot {tommy_lot} with number {tommy_number} "
    f"and trailer {tommy_trailer_number}. James owns Lot {james_lot} with "
    f"number {james_number} and trailer {james_trailer_number}."
)

# Although there is a legal issue with James's uncorrected deed,
# based on the recorded changes, this is the resulting ownership.
# Final Answer format
# print(f"<<<{final_answer}>>>")