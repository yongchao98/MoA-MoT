import os
import sys
import numpy as np
try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError:
    print("Pillow library not found. Please install it using: pip install Pillow")
    sys.exit(1)

def find_font(font_size=50):
    """
    Tries to find a font that supports Japanese characters on the system.
    """
    # Define a list of potential fonts and paths for different OSes.
    font_paths = [
        # Windows
        'C:/Windows/Fonts/meiryo.ttc',      # Meiryo
        'C:/Windows/Fonts/msgothic.ttc',    # MS Gothic
        'C:/Windows/Fonts/yugothr.ttc',     # Yu Gothic Regular
        # macOS
        '/System/Library/Fonts/Hiragino Sans GB.ttc', # Hiragino Sans
        '/System/Library/Fonts/Supplemental/Arial Unicode.ttf',
        # Linux (common locations)
        '/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',
        '/usr/share/fonts/truetype/takao-gothic/TakaoPGothic.ttf',
        '/usr/share/fonts/truetype/fonts-japanese-gothic.ttf',
        # Try loading by name, letting the system find it.
        'MS Gothic',
        'Meiryo',
        'Noto Sans CJK JP'
    ]
    for font_path in font_paths:
        try:
            return ImageFont.truetype(font_path, size=font_size)
        except (IOError, OSError):
            continue
    return None

def char_to_array(char, font, size=(64, 64)):
    """
    Renders a single character onto a fixed-size canvas and returns it as a numpy array.
    """
    # Create a black image
    image = Image.new('L', size, 0)
    draw = ImageDraw.Draw(image)

    # Calculate text bounding box and position to center it
    try:
        # Pillow >= 10.0.0 uses a 4-tuple bbox
        bbox = draw.textbbox((0, 0), char, font=font)
        text_width = bbox[2] - bbox[0]
        text_height = bbox[3] - bbox[1]
        # Position to draw the text to center its bounding box
        draw_pos = ((size[0] - text_width) / 2 - bbox[0], (size[1] - text_height) / 2 - bbox[1])
    except AttributeError:
        # Pillow < 10.0.0 uses textsize
        text_width, text_height = draw.textsize(char, font=font)
        draw_pos = ((size[0] - text_width) / 2, (size[1] - text_height) / 2)
    
    # Draw the character in white
    draw.text(draw_pos, char, font=font, fill=255)
    
    return np.array(image)

def jaccard_similarity(arr1, arr2):
    """
    Calculates the Jaccard similarity index between two boolean arrays.
    """
    # Convert arrays to boolean (True for any non-zero pixel)
    arr1_bool = arr1 > 128
    arr2_bool = arr2 > 128
    
    # Calculate intersection and union
    intersection = np.sum(np.logical_and(arr1_bool, arr2_bool))
    union = np.sum(np.logical_or(arr1_bool, arr2_bool))
    
    # Avoid division by zero if both images are blank
    if union == 0:
        return 1.0
    
    return intersection / union

def main():
    """
    Main function to run the visual comparison.
    """
    font = find_font()
    if not font:
        print("Error: Could not find a suitable Japanese font on your system.", file=sys.stderr)
        print("Please install a font like 'Noto Sans CJK' and try again.", file=sys.stderr)
        return

    target_char = 'ろ'
    candidate_digits = [str(i) for i in range(10)]

    print(f"Analyzing visual similarity for the Hiragana character: '{target_char}'")
    
    # Generate the image array for the target character
    target_array = char_to_array(target_char, font)
    
    similarities = {}
    
    # Compare the target character with each digit
    for digit in candidate_digits:
        digit_array = char_to_array(digit, font)
        similarity = jaccard_similarity(target_array, digit_array)
        similarities[digit] = similarity

    # Find the best match
    best_match_digit = max(similarities, key=similarities.get)

    print("\n--- Similarity Scores (higher is more similar) ---")
    # Sort by score for a clean report
    for digit, score in sorted(similarities.items(), key=lambda item: item[1], reverse=True):
        print(f"Similarity between 'ろ' and '{digit}': {score:.4f}")

    print("\n--- Conclusion ---")
    print(f"The Hiragana character 'ろ' is most visually similar to the digit: {best_match_digit}")

if __name__ == '__main__':
    main()