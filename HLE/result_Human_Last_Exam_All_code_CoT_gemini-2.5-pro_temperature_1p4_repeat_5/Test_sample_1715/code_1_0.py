# Initial information for Tommy and James
tommy_number = "1234567890"
tommy_trailer_number = "1234567890.1"
james_number = "0987654321"
james_trailer_number = "0987654321.1"

# Initial ownership in 2016
tommy_lot_2016 = "Lot A"
james_lot_2016 = "Lot B"

# In 2017, Tommy's Deed was legally corrected and re-recorded.
# This is the legally binding change.
tommy_lot_final = "Lot B"

# Since Tommy's Estate now legally owns Lot B, Lot A is left for James,
# regardless of whether he updated his own paperwork. The previous ownership
# state is legally invalidated by the corrected deed.
james_lot_final = "Lot A"


print("Final Ownership Status:")
print("=======================")
print(
    f"Tommy's Estate owns {tommy_lot_final} with Tax Map Number {tommy_number} "
    f"and trailer number {tommy_trailer_number}."
)
print(
    f"James owns {james_lot_final} with Tax Map Number {james_number} "
    f"and trailer number {james_trailer_number}."
)
