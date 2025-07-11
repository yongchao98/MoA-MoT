# This script identifies and displays the content of the Chinese seal.

def display_seal_content():
    """
    Identifies the characters on the seal, converts them to simplified Chinese,
    and prints the result with additional information.
    """
    # The seal is read from right to left, top to bottom.
    
    # Character identification based on position
    top_right = {
        "simplified": "量",
        "traditional": "量",
        "pinyin": "liáng"
    }
    bottom_right = {
        "simplified": "广",
        "traditional": "廣",
        "pinyin": "guǎng"
    }
    top_left = {
        "simplified": "天",
        "traditional": "天",
        "pinyin": "tiān"
    }
    bottom_left = {
        "simplified": "德",
        "traditional": "德",
        "pinyin": "dé"
    }

    # Assemble the phrase in simplified Chinese based on the reading order
    phrase = (top_right["simplified"] + 
              bottom_right["simplified"] + 
              top_left["simplified"] + 
              bottom_left["simplified"])
    
    print("The characters on the seal in simplified Chinese are:")
    print(f"Reading Order: Right-to-Left, Top-to-Bottom")
    print(f"Phrase: {phrase}")
    print(f"Pinyin: {top_right['pinyin']} {bottom_right['pinyin']} {top_left['pinyin']} {bottom_left['pinyin']}")
    print("\n--- Character Breakdown ---")
    print(f"1. Top Right: {top_right['simplified']} (pinyin: {top_right['pinyin']})")
    print(f"2. Bottom Right: {bottom_right['simplified']} (pinyin: {bottom_right['pinyin']}, traditional: {bottom_right['traditional']})")
    print(f"3. Top Left: {top_left['simplified']} (pinyin: {top_left['pinyin']})")
    print(f"4. Bottom Left: {bottom_left['simplified']} (pinyin: {bottom_left['pinyin']})")

if __name__ == '__main__':
    display_seal_content()