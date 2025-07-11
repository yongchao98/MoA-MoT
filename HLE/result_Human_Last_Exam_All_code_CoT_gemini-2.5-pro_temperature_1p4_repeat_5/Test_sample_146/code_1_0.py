import sys
from PIL import Image, ImageDraw, ImageFont
import numpy as np

def find_font():
    """
    Tries to find a suitable font on the system that can render Japanese characters.
    It checks common font names and paths for Windows, macOS, and Linux.
    """
    # A list of font filenames to try.
    font_names = [
        # Windows
        'meiryo.ttc', 'msmincho.ttc', 'msgothic.ttc', 'yumin.ttf',
        # macOS
        'Hiragino Sans GB W3.otf', 'Hiragino Sans W3.ttc',
        # Linux / Google Noto
        'NotoSansCJK-Regular.ttc', 'NotoSansJP-Regular.otf', 'DroidSansJapanese.ttf'
    ]
    # Common full paths
    font_paths = [
        '/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',
        '/System/Library/Fonts/Supplemental/Arial Unicode.ttf',
        'C:/Windows/Fonts/meiryo.ttc',
    ]

    for name in font_names:
        try:
            return ImageFont.truetype(name, 24)
        except IOError:
            continue
    for path in font_paths:
        try:
            return ImageFont.truetype(path, 24)
        except IOError:
            continue
            
    # As a last resort, try loading a default font. May not support 'ろ'.
    try:
        return ImageFont.load_default()
    except IOError:
        return None

def get_character_as_array(char, font, size=(32, 32)):
    """Renders a single character into a binary NumPy array."""
    # Create a white image
    image = Image.new('L', size, 255)
    draw = ImageDraw.Draw(image)
    
    # Center the character in the image
    try:
        # Modern Pillow versions use textbbox
        bbox = draw.textbbox((0, 0), char, font=font)
        text_width = bbox[2] - bbox[0]
        text_height = bbox[3] - bbox[1]
    except AttributeError:
        # Older Pillow versions use textsize
        text_width, text_height = draw.textsize(char, font=font)
        
    position = ((size[0] - text_width) / 2, (size[1] - text_height) / 2)
    
    # Draw the character in black
    draw.text(position, char, font=font, fill=0)
    
    # Convert to a NumPy array where True represents a black pixel
    return np.array(image) == 0

def calculate_jaccard_similarity(arr1, arr2):
    """Calculates the Jaccard similarity (intersection over union) for two boolean arrays."""
    intersection = np.sum(arr1 & arr2)
    union = np.sum(arr1 | arr2)
    return intersection / union if union > 0 else 0

def main():
    """
    Main function to find the most similar digit to 'ろ'.
    """
    # Ensure required libraries are installed
    try:
        from PIL import Image, ImageDraw, ImageFont
        import numpy as np
    except ImportError:
        print("Error: This script requires the 'Pillow' and 'numpy' libraries.")
        print("Please install them by running: pip install Pillow numpy")
        sys.exit(1)

    hiragana_char = 'ろ'
    digits = [str(i) for i in range(10)]

    # Find a working font
    font = find_font()
    if not font:
        print("Error: Could not find a suitable font to render Japanese characters.")
        print("Please ensure a font like 'Meiryo', 'Hiragino Sans', or 'Noto Sans CJK' is installed.")
        sys.exit(1)

    # Render the target character
    hiragana_array = get_character_as_array(hiragana_char, font)
    
    # Check if the font successfully rendered the character
    if np.sum(hiragana_array) == 0:
        print(f"Error: The located font failed to render the character '{hiragana_char}'.")
        print("Please install a Japanese-compatible font for this script to work.")
        sys.exit(1)

    best_match_digit = None
    max_similarity = -1.0
    
    print(f"Calculating visual similarity for hiragana '{hiragana_char}' against digits 0-9.")
    print("--------------------------------------------------")

    # Loop through each digit, render it, and compare it to 'ろ'
    for digit in digits:
        digit_array = get_character_as_array(digit, font)
        similarity = calculate_jaccard_similarity(hiragana_array, digit_array)
        
        # This is the "equation" outputting the similarity score for each digit
        print(f"Similarity equation: Jaccard('{hiragana_char}', '{digit}') = {similarity:.4f}")
        
        if similarity > max_similarity:
            max_similarity = similarity
            best_match_digit = digit
    
    print("--------------------------------------------------")
    print(f"The digit most visually similar to '{hiragana_char}' is '{best_match_digit}'.")

if __name__ == "__main__":
    main()