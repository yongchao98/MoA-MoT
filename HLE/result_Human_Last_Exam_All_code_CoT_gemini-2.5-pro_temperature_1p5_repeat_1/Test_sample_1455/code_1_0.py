import math

def solve_counterpoint():
    """
    Analyzes the West Side Story counterpoint excerpt to find X, Y, and Z.
    """
    # Step 1 & 2: Define and represent the melodies with MIDI numbers.
    # REST is represented by 0. C4 is 60.
    # The melodic lines have been determined based on the clues and principle of minimal melodic leaps.
    # A key finding was that 3-event measures follow a 3+2+1 beat rhythm.
    maria_melody = [
        # M0: x x x E E E (1+1+1+1+1+1)
        0, 0, 0, 64, 64, 64,
        # M1: A---- G-- A (3+2+1)
        69, 69, 69, 67, 67, 69,
        # M2: B♭---- C---- (3+3)
        70, 70, 70, 72, 72, 72,
        # M3: F---- C---- (3+3)
        65, 65, 65, 72, 72, 72,
        # M4: F---- F---- (3+3)
        65, 65, 65, 65, 65, 65,
        # M5: G---- A---- (3+3)
        67, 67, 67, 69, 69, 69,
        # M6: B♮---- C----- (3+3)
        71, 71, 71, 72, 72, 72,
        # M7: ------------ (6)
        0, 0, 0, 0, 0, 0
    ]

    tony_melody = [
        # M0: x x x E G C (1+1+1+1+1+1)
        0, 0, 0, 64, 55, 60,
        # M1: F------ G F (3+2+1)
        65, 65, 65, 67, 67, 65,
        # M2: D---- E---- (3+3)
        62, 62, 62, 64, 64, 64,
        # M3: D---- A---- (3+3)
        62, 62, 62, 69, 69, 69,
        # M4: D---- F-- G (3+2+1)
        62, 62, 62, 65, 65, 67,
        # M5: E-- F D----- (3+2+1)
        64, 64, 64, 65, 65, 62,
        # M6: --- E♭ C----- (3+2+1)
        0, 0, 0, 63, 63, 60,
        # M7: ------------ (6)
        0, 0, 0, 0, 0, 0
    ]

    # Step 3.1: Find X and Y
    # Find the single beat where Tony's voice is higher than Maria's.
    x_val, y_val = -1, -1
    higher_beats_count = 0
    for i in range(len(maria_melody)):
        # Both must be singing
        if tony_melody[i] > 0 and maria_melody[i] > 0:
            if tony_melody[i] > maria_melody[i]:
                higher_beats_count += 1
                measure = i // 6
                beat = i % 6
                x_val, y_val = measure, beat

    # Step 3.2: Find Z
    # Count instances of contrary motion.
    z_val = 0
    for i in range(1, len(maria_melody)):
        m_curr, m_prev = maria_melody[i], maria_melody[i-1]
        t_curr, t_prev = tony_melody[i], tony_melody[i-1]

        # Condition: Both must be singing on the current and previous beat.
        if all(n > 0 for n in [m_curr, m_prev, t_curr, t_prev]):
            # Condition: Both melodies must move.
            if m_curr != m_prev and t_curr != t_prev:
                # Get direction of motion (-1 for down, 1 for up)
                maria_motion = math.copysign(1, m_curr - m_prev)
                tony_motion = math.copysign(1, t_curr - t_prev)
                # Check if directions are opposite
                if maria_motion * tony_motion == -1:
                    z_val += 1

    print(f"{x_val} {y_val} {z_val}")

solve_counterpoint()