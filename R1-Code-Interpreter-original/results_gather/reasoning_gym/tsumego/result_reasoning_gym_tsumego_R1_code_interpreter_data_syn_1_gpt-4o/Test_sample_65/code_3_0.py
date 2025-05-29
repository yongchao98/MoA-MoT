# Simulate the move sequence to verify the capture
def simulate_capture():
    # Initial liberties for the white stones at D6 and E8
    liberties_D6 = 3  # C6, D5, and D7
    liberties_E8 = 3  # D7, E7, and F8

    # Black plays at D8
    liberties_D6 -= 1  # D8 is occupied
    liberties_E8 -= 1  # D8 is occupied

    # Check if the stones at D6 and E8 are captured
    if liberties_D6 == 0 or liberties_E8 == 0:
        return "Capture successful at D8"
    else:
        return "Capture failed"

print(simulate_capture())