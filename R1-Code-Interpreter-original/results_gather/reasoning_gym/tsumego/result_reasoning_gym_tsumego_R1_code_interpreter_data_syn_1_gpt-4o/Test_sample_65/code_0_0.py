# Simulate the move sequence to verify the capture
def simulate_capture():
    # Initial liberties for the white stones at D6 and E8
    liberties_D6 = 2  # C6 and D5
    liberties_E8 = 3  # D7, E7, and F8

    # Black plays at C6
    liberties_D6 -= 1  # C6 is occupied
    if liberties_D6 == 1:
        # White plays at C5 to escape
        liberties_D6 += 1  # C5 is now a liberty

    # Black plays at D5 to capture
    liberties_D6 -= 1  # D5 is occupied

    # Check if the stones at D6 and E8 are captured
    if liberties_D6 == 0:
        return "Capture successful at C6"
    else:
        return "Capture failed"

print(simulate_capture())