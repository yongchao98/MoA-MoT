import math

def solve_west_side_story():
    """
    This function solves the West Side Story counterpoint puzzle by programmatically
    analyzing the derived pitch data for the two vocal lines.
    """
    
    # Step 1 & 2: The musical notation is parsed and pitches are assigned based on the clues.
    # The following pitch assignments (in MIDI numbers) are the unique solution that
    # satisfies all the given constraints.
    # None represents a rest.
    
    maria_pitches = [
        None, None, None, 64, 64, 64,  # M0: x x x E E E (E4)
        69, 69, 69, 69, 67, 67,          # M1: A---- G-- (A4, G4)
        70, 70, 70, 70, 72, 72,          # M2: B♭---- C-- (Bb4, C5)
        65, 65, 65, 65, 60, 60,          # M3: F---- C-- (F4, C4)
        65, 65, 65, 65, 65, 65,          # M4: F---- F---- (F4)
        67, 67, 67, 67, 69, 69,          # M5: G---- A-- (G4, A4)
        71, 71, 71, 71, 72, 72,          # M6: B♮---- C-- (B4, C5)
        None, None, None, None, None, None # M7: rest
    ]

    tony_pitches = [
        None, None, None, 64, 55, 60,  # M0: x x x E G C (E4, G3, C4)
        65, 65, 65, 65, 65, 65,          # M1: F------ (F4)
        62, 62, 62, 62, 64, 64,          # M2: D---- E-- (D4, E4)
        62, 62, 62, 62, 57, 57,          # M3: D---- A-- (D4, A3)
        62, 62, 62, 62, 65, 67,          # M4: D---- F G (D4, F4, G4)
        64, 64, 65, 65, 62, 62,          # M5: E-- F-- D-- (E4, F4, D4)
        None, None, None, 63, 63, 60,  # M6: --- E♭-- C (rest, Eb4, C4)
        None, None, None, None, None, None # M7: rest
    ]

    # Step 3: Determine X (measure) and Y (beat)
    X, Y = -1, -1
    for i in range(len(maria_pitches)):
        # Check for beats where both are singing
        if maria_pitches[i] is not None and tony_pitches[i] is not None:
            if tony_pitches[i] > maria_pitches[i]:
                X = i // 6
                Y = i % 6
                # Since the prompt says this happens for a total of one beat, we can stop.
                break

    # Step 4: Determine Z (contrary motion count)
    Z = 0
    # First, create a list of "events" which are points where the pitches change.
    # An event is defined by the pitches of both singers.
    events = []
    last_m_pitch = None
    last_t_pitch = None
    for i in range(len(maria_pitches)):
        m_pitch = maria_pitches[i]
        t_pitch = tony_pitches[i]
        
        # An event is a change in pitch for either singer.
        is_new_event = (m_pitch != last_m_pitch or t_pitch != last_t_pitch)
        
        # We only count motion when both singers are actively singing at the event point.
        if is_new_event and m_pitch is not None and t_pitch is not None:
            events.append({'m': m_pitch, 't': t_pitch})
            
        last_m_pitch = m_pitch
        last_t_pitch = t_pitch

    # Compare consecutive events to find contrary motion
    for i in range(1, len(events)):
        m_prev = events[i-1]['m']
        t_prev = events[i-1]['t']
        m_curr = events[i]['m']
        t_curr = events[i]['t']

        m_delta = m_curr - m_prev
        t_delta = t_curr - t_prev

        # Contrary motion is when one voice ascends and the other descends.
        # This is true if the product of their pitch deltas is negative.
        if m_delta * t_delta < 0:
            Z += 1

    # Step 5: Print the final answer
    print(f"{X} {Y} {Z}")

solve_west_side_story()
<<<4 5 3>>>