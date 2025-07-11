import sys
import numpy as np
from PIL import Image, ImageDraw, ImageFont

def render_char_to_array(char_to_render, font, size=(64, 64)):
    """
    Renders a single character into a binary NumPy array.
    """
    # Create a blank, grayscale image
    image = Image.new("L", size, "white")
    draw = ImageDraw.Draw(image)
    
    # Calculate text size to center it on the canvas
    try:
        # Newer Pillow versions use textbbox
        bbox = draw.textbbox((0, 0), char_to_render, font=font)
        text_width = bbox[2] - bbox[0]
        text_height = bbox[3] - bbox[1]
    except AttributeError:
        # Older Pillow versions use textsize
        text_width, text_height = draw.textsize(char_to_render, font=font)
    
    position = ((size[0] - text_width) / 2, (size[1] - text_height) / 2)
    
    # Draw the character in black
    draw.text(position, char_to_render, font=font, fill="black")
    
    # Convert the image to a NumPy array and binarize it
    # 1 represents a black pixel (the character), 0 represents a white pixel (background)
    np_image = np.array(image)
    return (np_image < 128).astype(int)

def find_font(font_size=56):
    """
    Tries to find a suitable font on the system that supports Japanese.
    """
    # Common font paths for different operating systems
    # We prioritize Japanese fonts for the best rendering of 'ろ'.
    font_paths = [
        # Windows
        'C:/Windows/Fonts/yugothb.ttc',   # Yu Gothic
        'C:/Windows/Fonts/meiryo.ttc',    # Meiryo
        'C:/Windows/Fonts/msgothic.ttc',  # MS Gothic
        # macOS
        '/System/Library/Fonts/ヒラギノ角ゴ ProN.ttc', # Hiragino Sans
        '/System/Library/Fonts/Supplemental/Arial Unicode.ttf',
        # Linux
        '/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc', # Noto Sans CJK
        '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf' # Fallback
    ]
    
    for path in font_paths:
        try:
            return ImageFont.truetype(path, font_size)
        except IOError:
            continue
            
    # If no preferred font is found, try to load a default
    try:
        return ImageFont.load_default()
    except IOError:
        return None

def main():
    """
    Main function to run the character similarity analysis.
    """
    font = find_font()
    if font is None:
        print("Error: Could not find a suitable font on your system.")
        print("Please ensure you have a standard CJK or Unicode font installed.")
        sys.exit(1)

    target_char = 'ろ'
    digits = '0123456789'
    
    # Render the target character 'ろ'
    ro_array = render_char_to_array(target_char, font)
    
    similarities = {}
    
    print(f"Analyzing which digit is most visually similar to '{target_char}':\n")

    # Compare 'ろ' with each digit
    for digit in digits:
        digit_array = render_char_to_array(digit, font)
        
        # Calculate Sum of Absolute Differences (SAD)
        # This counts the number of pixels that are different
        difference = np.sum(np.abs(ro_array - digit_array))
        similarities[digit] = difference
        
        # Print the "equation" showing the comparison for each number
        print(f"Similarity of '{target_char}' vs '{digit}': Difference Score = {difference}")

    # Find the digit with the minimum difference score
    most_similar_digit = min(similarities, key=similarities.get)
    min_difference = similarities[most_similar_digit]
    
    print("\n---")
    print("Conclusion:")
    print(f"The character '{most_similar_digit}' has the lowest difference score ({min_difference}).")
    print(f"Therefore, the hiragana character 'ろ' is most likely to be mistaken for the digit {most_similar_digit}.")

if __name__ == '__main__':
    main()