import math

# Step 1: Determine the value of S
S_str = "SNYWR"
BC = len(S_str)  # Side a, adjacent to angle B

# Step 2: Define the given angle
B_deg = 30  # Angle B in degrees

# Step 3: Explain the geometric principle (Ambiguous Case of Law of Sines)
print(f"The value of S = len('SNYWR') is {BC}.")
print(f"In triangle ABC, side BC = {BC} and angle B = {B_deg} degrees.")
print("For there to be two possible values for side AB, we must be in the ambiguous case (SSA) of the Law of Sines.")
print("This occurs when the side opposite the angle (AC) is greater than the altitude from vertex C (h) but less than the other given side (BC).")

# Step 4: Calculate the altitude h
B_rad = math.radians(B_deg)
h = BC * math.sin(B_rad)
print(f"\nThe altitude h = BC * sin(B) = {BC} * sin({B_deg}) = {h:.1f}")

# Step 5: Establish and print the inequality for two triangles
print("\nThe condition for two triangles is: h < AC < BC")
print("Substituting the values, the final inequality is:")
print(f"{h:.1f} < AC < {BC}")

# Step 6: Find and count the integer values for AC that satisfy the inequality
# The smallest integer greater than h
start_int = math.floor(h) + 1
# The largest integer less than BC
end_int = BC - 1

integer_values = list(range(start_int, end_int + 1))
count = len(integer_values)

print(f"\nThe integer values for AC that satisfy this condition are: {', '.join(map(str, integer_values))}.")
print(f"\nThe number of integer values is {count}.")