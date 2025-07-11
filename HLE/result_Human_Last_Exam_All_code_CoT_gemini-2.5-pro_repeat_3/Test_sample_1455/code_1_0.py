import re

def parse_line(line, octave_map):
    """Parses a measure string into a list of 6 pitches."""
    # This parser assumes note names are 1 or 2 chars (e.g., C, B♭)
    # and are followed by dashes for duration.
    # e.g., "A---- G--" -> A for 4 beats, G for 2 beats.
    # e.g., "x x x E G C" -> rest, rest, rest, E, G, C
    
    note_to_midi_base = {
        'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11
    }
    
    beats = [None] * 6
    current_beat = 0
    
    parts = line.split()
    for part in parts:
        if part == 'x':
            current_beat += 1
            continue
        
        # Extract note name and accidental
        note_name = part[0]
        accidental = ''
        if len(part) > 1 and part[1] in '♭♮':
            accidental = part[1]
            note_key = note_name + accidental
        else:
            note_key = note_name

        # Determine duration
        duration = len(part.replace(note_key, '', 1)) + 1
        
        # Get base pitch and apply octave
        base_pitch = note_to_midi_base[note_name]
        if accidental == '♭':
            base_pitch -= 1
        
        octave = octave_map.get(note_key, 4) # Default to octave 4 if not specified
        pitch = base_pitch + 12 * (octave + 1)
        
        for i in range(duration):
            if current_beat < 6:
                beats[current_beat] = pitch
                current_beat += 1
                
    return beats

def solve_counterpoint():
    """
    Solves the West Side Story counterpoint puzzle.
    """
    maria_raw = [
        "x x x E E E", "A---- G--", "B♭---- C----", "F---- C----",
        "F---- F----", "G---- A----", "B♮---- C-----", "x x x x x x"
    ]
    tony_raw = [
        "x x x E G C", "F---- G F", "D---- E----", "D---- A----",
        "D---- F-- G", "E-- F D-----", "--- E♭ C-----", "x x x x x x"
    ]
    
    # Octave assignments deduced from the clues
    # Format: {NoteName: OctaveNumber}
    maria_octaves = {'C': 5, 'B♮': 4, 'B♭': 4, 'A': 4, 'G': 4, 'F': 4, 'E': 4}
    
    # This specific assignment for Tony satisfies all clues, including "lowest note is G"
    tony_octaves = {
        'C': 4, 'A': 4, 'G': 4, 'F': 4, 'E': 4, 'E♭': 4, 'D': 4,
        # Special cases to satisfy all constraints
        ('G', 0): 3, # G in measure 0 is G3
    }
    
    # Manually build the pitch arrays based on deductions
    # This ensures perfect accuracy according to the problem's constraints
    m_pitches = [
        None, None, None, 64, 64, 64,  # M0
        69, 69, 69, 69, 67, 67,      # M1
        70, 70, 70, 70, 72, 72,      # M2
        65, 65, 65, 65, 72, 72,      # M3
        65, 65, 65, 65, 65, 65,      # M4
        67, 67, 67, 67, 69, 69,      # M5
        71, 71, 71, 71, 72, 72,      # M6
        None, None, None, None, None, None # M7
    ]
    
    t_pitches = [
        None, None, None, 64, 55, 60,  # T0 (E4, G3, C4)
        65, 65, 65, 65, 67, 65,      # T1
        62, 62, 62, 62, 64, 64,      # T2
        62, 62, 62, 62, 69, 69,      # T3
        62, 62, 62, 62, 65, 67,      # T4
        64, 64, 65, 62, 62, 62,      # T5
        None, None, None, 63, 60, 60,      # T6
        None, None, None, None, None, None # T7
    ]
    
    # --- Find X and Y ---
    X, Y = -1, -1
    for i in range(len(m_pitches)):
        if m_pitches[i] is not None and t_pitches[i] is not None:
            if t_pitches[i] > m_pitches[i]:
                X = i // 6
                Y = i % 6
                break # Found the single occurrence

    # --- Find Z ---
    Z = 0
    # Get a list of times where at least one voice sings a new note
    event_times = set()
    for i in range(1, len(m_pitches)):
        if m_pitches[i] != m_pitches[i-1] or t_pitches[i] != t_pitches[i-1]:
            event_times.add(i)
    
    sorted_times = sorted(list(event_times))

    # To handle the first note
    if m_pitches[0] is not None or t_pitches[0] is not None:
         sorted_times.insert(0,0)
    
    # Check motion between events
    for i in range(len(sorted_times) - 1):
        start_time = sorted_times[i]
        end_time = sorted_times[i+1]
        
        m_start, m_end = m_pitches[start_time], m_pitches[end_time]
        t_start, t_end = t_pitches[start_time], t_pitches[end_time]

        # Both voices must be singing at both ends of the interval
        if all(p is not None for p in [m_start, m_end, t_start, t_end]):
            dM = m_end - m_start
            dT = t_end - t_start
            
            # Both voices must move
            if dM != 0 and dT != 0:
                # Check for contrary motion
                if (dM > 0 and dT < 0) or (dM < 0 and dT > 0):
                    Z += 1

    print(f"{X} {Y} {Z}")

solve_counterpoint()