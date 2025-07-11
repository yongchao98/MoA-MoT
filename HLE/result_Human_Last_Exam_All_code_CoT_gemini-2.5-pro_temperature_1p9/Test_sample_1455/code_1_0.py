import re

def solve_counterpoint():
    """
    Solves the West Side Story counterpoint puzzle.
    """
    maria_raw = "x x x E E E |A---- G-- A |B♭---- C---- |F---- C---- |F---- F---- |G---- A---- |B♮---- C-----|------------"
    tony_raw = "x x x E G C |F------ G F |D---- E---- |D---- A---- |D---- F-- G |E-- F D-----|--- E♭ C-----|------------"

    NOTE_TO_BASE_MIDI = {
        'C': 0, 'B♮': 11, 'B': 11,
        'B♭': 10, 'A': 9, 'G': 7,
        'F': 5, 'E': 4, 'E♭': 3,
        'D': 2, 'x': -1, '-': -1
    }

    def parse_note_token(token):
        """Extracts note name from a token like 'B♭----' -> 'B♭'"""
        if token.startswith('---'): return '-'
        if token.startswith('x'): return 'x'
        match = re.match(r"([A-G][b♮♭]?)", token)
        return match.group(1) if match else None

    def get_beat_notes(raw_line):
        """Parses the raw string into a list of 48 beat notes."""
        measures_str = [s.strip() for s in raw_line.split('|')]
        beat_notes = []
        for measure_str in measures_str:
            tokens = measure_str.split()
            notes_in_measure = []
            if not tokens: # Handle empty measure like | |
                notes_in_measure = ['-'] * 6
            else:
                num_tokens = len(tokens)
                if num_tokens == 6: # 1 beat each
                    durations = [1] * 6
                elif num_tokens == 3: # 2 beats each
                    durations = [2] * 3
                elif num_tokens == 2: # 3 beats each
                    durations = [3] * 2
                elif num_tokens == 1: # 6 beats
                    durations = [6]
                else: # Fallback for unexpected token counts
                    durations = [6 // num_tokens] * num_tokens
                
                for i, token in enumerate(tokens):
                    note = parse_note_token(token)
                    notes_in_measure.extend([note] * durations[i])
            beat_notes.extend(notes_in_measure[:6]) # Ensure exactly 6 beats
        return beat_notes

    def resolve_pitches(beat_notes):
        """Converts note names to MIDI pitches, resolving octaves."""
        pitches = []
        last_pitch = 0
        
        # Start unison 'E' on E4, which fits all constraints except one.
        initial_pitch_context = 64 

        for note in beat_notes:
            if note is None or note in ('x', '-'):
                pitches.append(0)
                continue

            base_midi = NOTE_TO_BASE_MIDI[note]
            
            # First note
            if last_pitch == 0:
                # Set initial E to E4 (64)
                current_pitch = initial_pitch_context if note == 'E' else None
                if current_pitch is None: # Should not happen after first note
                     # Failsafe: find a reasonable starting octave
                     current_pitch = 50 + base_midi
            else:
                # Find the closest octave for the new note
                # Check options one octave down, same octave, one octave up
                options = []
                last_octave = last_pitch - (last_pitch % 12)
                for octave_offset in [-12, 0, 12]:
                    candidate = last_octave + octave_offset + base_midi
                    # Also consider case where C is above B
                    if last_pitch % 12 > base_midi and octave_offset == 0:
                         candidate += 12
                    # Consider case where B is below C
                    if last_pitch % 12 < base_midi and octave_offset == 0:
                        candidate -=12 if (last_octave + base_midi > last_pitch + 9) else 0

                    options.append(candidate)
                
                # Choose the option with the smallest jump (<= Major 6th / 9 semitones)
                valid_options = [p for p in options if abs(p - last_pitch) <= 9]
                if not valid_options:
                     # If no jump fits, just take the closest
                     valid_options = options
                
                current_pitch = min(valid_options, key=lambda p: abs(p - last_pitch))

            pitches.append(current_pitch)
            if current_pitch != 0:
                last_pitch = current_pitch
        return pitches
    
    # Run the parsing and analysis
    m_notes = get_beat_notes(maria_raw)
    t_notes = get_beat_notes(tony_raw)
    
    # Combine notes and resolve octaves in tandem to maintain relative positions
    combined_notes = list(zip(m_notes, t_notes))
    m_pitches, t_pitches = [], []
    last_m_pitch, last_t_pitch = 0, 0
    # Special handling for unison start
    m_pitches = [0,0,0,64,64,64]
    t_pitches = [0,0,0,64,67,72]
    last_m_pitch, last_t_pitch = 64, 72

    # A more deterministic octave resolver based on the analysis
    def get_pitch(note, prev_pitch):
        if not note or note in ('x', '-'): return 0
        base = NOTE_TO_BASE_MIDI[note]
        # Check candidate octaves
        candidates = []
        for i in range(-2, 3):
            candidates.append(prev_pitch - (prev_pitch % 12) + i * 12 + base)
        
        # Adjust for wrapping around C
        if base < prev_pitch % 12:
            candidates.append(prev_pitch - (prev_pitch % 12) + 12 + base)
        else:
            candidates.append(prev_pitch - (prev_pitch % 12) -12 + base)
        
        valid_candidates = [c for c in candidates if abs(c - prev_pitch) <= 9] if prev_pitch != 0 else candidates
        if not valid_candidates: valid_candidates = candidates # Relax if jump is too big
        
        return min(valid_candidates, key=lambda c: abs(c - prev_pitch))

    for i in range(6, 48):
        m_note, t_note = combined_notes[i]

        m_pitch = get_pitch(m_note, last_m_pitch)
        t_pitch = get_pitch(t_note, last_t_pitch)
        
        # Manual check to ensure they are within an octave
        if m_pitch != 0 and t_pitch != 0 and abs(m_pitch-t_pitch) > 12:
            # Re-calculate Tony's note relative to Maria's
            t_pitch_alt = get_pitch(t_note, m_pitch)
            if abs(t_pitch_alt - last_t_pitch) <= 9:
                t_pitch = t_pitch_alt

        m_pitches.append(m_pitch)
        t_pitches.append(t_pitch)
        if m_pitch != 0: last_m_pitch = m_pitch
        if t_pitch != 0: last_t_pitch = t_pitch

    # Find X and Y
    X, Y = -1, -1
    for i in range(48):
        if t_pitches[i] > m_pitches[i]:
            X = i // 6
            Y = i % 6
            break
            
    # Find Z
    Z = 0
    m_last_note, t_last_note = 0, 0
    for i in range(48):
        m_curr, t_curr = m_pitches[i], t_pitches[i]
        
        m_moved = (m_curr != m_last_note and m_curr != 0 and m_last_note != 0)
        t_moved = (t_curr != t_last_note and t_curr != 0 and t_last_note != 0)
        
        if m_moved and t_moved:
            m_dir = 1 if m_curr > m_last_note else -1
            t_dir = 1 if t_curr > t_last_note else -1
            if m_dir * t_dir == -1:
                Z += 1

        if m_curr != 0: m_last_note = m_curr
        if t_curr != 0: t_last_note = t_curr
        
    print(f"{X} {Y} {Z}")

solve_counterpoint()