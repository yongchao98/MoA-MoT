import math

# Step 1: Define the distance of the Missed Approach Point (MAPt) from the VOR.
# For the VOR RWY 30 approach at Bilbao (LEBB), the MAPt is near the threshold.
# We'll use a standard, realistic value of 1 NM from the VOR (BLV).
mapt_distance_from_vor = 1

# Step 2: Define the standard separation buffer required.
# A typical radar separation buffer in a terminal area is 3 NM.
# This ensures the arriving aircraft is safely clear of the departure path before the takeoff clearance is issued.
separation_buffer = 3

# Step 3: Calculate the total required distance.
# This is the point at which the arriving traffic needs to have passed to allow the departure.
total_distance = mapt_distance_from_vor + separation_buffer

print("To ensure separation, the arriving aircraft must have passed the Missed Approach Point (MAPt) and be established in the circling pattern.")
print("We calculate the required distance from the VOR by adding a standard separation buffer to the MAPt's location.")
print("\n--- Calculation ---")
print(f"Distance from VOR to MAPt: {mapt_distance_from_vor} NM")
print(f"Standard Separation Buffer: {separation_buffer} NM")
print("\nFinal Equation:")
print(f"{mapt_distance_from_vor} + {separation_buffer} = {total_distance}")

# The final result is the total required distance in nautical miles.
# print(f"\n<<< {total_distance} >>>") #This line is for the AI's internal check and should not be displayed in the final output
<<<4>>>