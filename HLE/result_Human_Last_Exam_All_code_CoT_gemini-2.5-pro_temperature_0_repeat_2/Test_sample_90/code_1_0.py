import librosa
import numpy as np

def transcribe_melody(file_path, start_sec, end_sec):
    """
    Analyzes a segment of an audio file to identify the melody notes.

    Args:
        file_path (str): The path to the audio file (e.g., .mp3, .wav).
        start_sec (int): The start time of the segment in seconds.
        end_sec (int): The end time of the segment in seconds.
    """
    try:
        # Load the specified segment of the audio file
        # sr=None preserves the original sample rate
        y, sr = librosa.load(file_path, offset=start_sec, duration=end_sec - start_sec, sr=None)

        # Use the PYIN algorithm to estimate the fundamental frequency (pitch)
        # This is effective for melodic instruments like a piano's right hand.
        # We set a frequency range typical for a piano melody.
        f0, voiced_flag, voiced_probs = librosa.pyin(
            y,
            fmin=librosa.note_to_hz('C3'),
            fmax=librosa.note_to_hz('C6')
        )

        # Filter out timestamps where no pitch was detected (NaN values)
        pitches_hz = f0[~np.isnan(f0)]

        if len(pitches_hz) == 0:
            print("Could not detect any notes in the specified time frame.")
            print("Please check the file path and the time segment.")
            return

        # Convert the frequencies (in Hz) to musical note names
        # e.g., 440.0 -> 'A4'
        detected_notes = librosa.hz_to_note(pitches_hz, sharps=True)

        # The pitch detection creates many notes. We'll simplify the list
        # to show the melodic progression by removing consecutive duplicates.
        melody_line = []
        if len(detected_notes) > 0:
            # Add the first note
            melody_line.append(detected_notes[0])
            # Add subsequent notes only if they are different from the previous one
            for i in range(1, len(detected_notes)):
                if detected_notes[i] != detected_notes[i-1]:
                    melody_line.append(detected_notes[i])
        
        # Format the output to use ♭ for flats instead of the default # for sharps
        # This is more common in some musical contexts.
        note_map = {
            'A#': 'B♭', 'C#': 'D♭', 'D#': 'E♭', 'F#': 'G♭', 'G#': 'A♭'
        }
        formatted_melody = []
        for note in melody_line:
            pitch_class = note[:-1]
            octave = note[-1]
            if pitch_class in note_map:
                formatted_melody.append(note_map[pitch_class] + octave)
            else:
                formatted_melody.append(note)


        print("Detected melody notes from {}s to {}s:".format(start_sec, end_sec))
        # Print each note in the final sequence
        print(' '.join(formatted_melody))
        print("\nNote: This is an automated transcription and may not be 100% accurate.")

    except FileNotFoundError:
        print(f"Error: The file was not found at '{file_path}'")
        print("Please make sure the file path is correct.")
    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have 'ffmpeg' installed and accessible in your system's PATH.")
        print("You may also need to install libraries with: pip install librosa numpy")


if __name__ == '__main__':
    # --- USER INPUT ---
    # IMPORTANT: Replace this with the actual path to your audio file.
    audio_file = 'path/to/your/song.mp3'
    
    # Time segment from the song "All My Loves Are You"
    start_time = 39
    end_time = 45

    transcribe_melody(audio_file, start_time, end_time)