# Simulate the move sequence to verify the capture
def simulate_capture():
    # Initial liberties for the white stones at D6 and E8
    liberties_D6 = 3  # C6, D5, and E7
    liberties_E8 = 3  # D7, E7, and F8

    # Black plays at E7
    liberties_D6 -= 1  # E7 is occupied
    liberties_E8 -= 1  # E7 is occupied

    # Check if the stones at D6 and E8 are captured
    if liberties_D6 == 0 or liberties_E8 == 0:
        return "Capture successful at E7"
    else:
        return "Capture failed"

print(simulate_capture())