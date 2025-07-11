import math

# The fraction we want to represent in Red-Blue Hackenbush
numerator = 13
denominator = 16

print(f"Finding the number of pieces for the fraction {numerator}/{denominator} in Red-Blue Hackenbush.\n")
print("The process involves constructing a 'stalk' of colored pieces recursively.")
print("At each step, we check if the current value is positive or negative to determine the piece's color, then calculate the next value.\n")

# Make copies of the numerator and denominator to modify in the loop
current_num = numerator
current_den = denominator

# Initialize a counter for the number of pieces
piece_count = 0

# Loop until the value of the fraction becomes 0
while current_num != 0:
    piece_count += 1
    
    # Store the fraction before updating for printing
    prev_num, prev_den = current_num, current_den

    # Determine the color and the next value based on the sign of the current value
    if current_num > 0:
        color = "Blue"
        # The new value is 2x - 1, which corresponds to a new numerator of 2*num - den
        next_num = 2 * current_num - current_den
        next_den = current_den
        
        print(f"Step {piece_count}:")
        print(f"  Value is {prev_num}/{prev_den} (Positive). The piece is {color}.")
        print(f"  The new value is calculated as 2 * ({prev_num}/{prev_den}) - 1.")
        print(f"  Equation: (2 * {prev_num} - {prev_den}) / {prev_den} = {next_num}/{next_den}")

    else: # This means current_num < 0
        color = "Red"
        # The new value is 2x + 1, which corresponds to a new numerator of 2*num + den
        next_num = 2 * current_num + current_den
        next_den = current_den

        print(f"Step {piece_count}:")
        print(f"  Value is {prev_num}/{prev_den} (Negative). The piece is {color}.")
        print(f"  The new value is calculated as 2 * ({prev_num}/{prev_den}) + 1.")
        print(f"  Equation: (2 * {prev_num} + {prev_den}) / {prev_den} = {next_num}/{next_den}")
    
    # Update the numerator and denominator for the next iteration
    current_num = next_num
    current_den = next_den

    # For display purposes, we can simplify the fraction for the next step
    if current_num != 0:
        common_divisor = math.gcd(current_num, current_den)
        current_num //= common_divisor
        current_den //= common_divisor
    
    print("-" * 30)

print(f"\nThe process stops because the value is now 0.")
print(f"The total number of pieces needed is {piece_count}.")

print(f"<<<{piece_count}>>>")