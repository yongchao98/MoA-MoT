import numpy as np
from PIL import Image, ImageDraw, ImageFont
from skimage.metrics import structural_similarity as ssim
import os

def char_to_image(char, font, size=(64, 64), mode='L'):
    """
    Renders a single character to a grayscale image and returns it as a NumPy array.
    """
    # Create a blank image
    image = Image.new(mode, size, color='white')
    draw = ImageDraw.Draw(image)

    # Get text size
    try:
        # Use getbbox for modern Pillow versions
        bbox = draw.textbbox((0, 0), char, font=font)
        text_width = bbox[2] - bbox[0]
        text_height = bbox[3] - bbox[1]
    except AttributeError:
        # Fallback for older Pillow versions
        text_width, text_height = draw.textsize(char, font=font)

    # Position the text
    x = (size[0] - text_width) / 2
    y = (size[1] - text_height) / 2
    
    # Draw the character onto the image
    draw.text((x, y), char, font=font, fill='black')
    
    # Convert image to a numpy array for comparison
    return np.array(image)

def find_most_similar_digit():
    """
    Finds which digit (0-9) is most visually similar to the hiragana 'ろ'.
    """
    # --- IMPORTANT ---
    # Please replace this with the path to a Japanese font file on your system.
    # For example, on Windows, it might be 'C:/Windows/Fonts/meiryo.ttc'
    # On macOS, it might be '/System/Library/Fonts/Hiragino Sans GB.ttc'
    # If you don't have one, you can download a free one like "Noto Sans JP".
    font_path = None #<-- SET YOUR FONT PATH HERE e.g., "C:/Windows/Fonts/YuGothM.ttc"
    
    # List of common font file names for different OS
    font_files = {
        'win32': ['meiryo.ttc', 'yugothm.ttc', 'msgothic.ttc'],
        'darwin': ['Hiragino Sans GB.ttc', 'Hiragino Mincho ProN.ttc'],
        'linux': ['NotoSansCJK-Regular.ttc', 'DroidSansJapanese.ttf']
    }
    
    # Common font directories
    font_dirs = {
        'win32': 'C:/Windows/Fonts',
        'darwin': '/System/Library/Fonts',
        'linux': ['/usr/share/fonts/truetype', '/usr/share/fonts/opentype']
    }
    
    import sys
    if sys.platform in font_dirs:
        dirs_to_check = font_dirs[sys.platform]
        if not isinstance(dirs_to_check, list):
            dirs_to_check = [dirs_to_check]
        
        for d in dirs_to_check:
            if not font_path and os.path.isdir(d):
                for font_name in font_files[sys.platform]:
                    potential_path = os.path.join(d, font_name)
                    if os.path.exists(potential_path):
                        font_path = potential_path
                        break
            if font_path:
                break
    
    if not font_path:
        print("Warning: A suitable Japanese font was not automatically found.")
        print("Please set the 'font_path' variable in the script manually.")
        # As a last resort, try loading a common name and hope PIL finds it
        font_path = 'Arial Unicode MS'


    try:
        font = ImageFont.truetype(font_path, 48)
    except IOError:
        print(f"Error: The font file could not be found or loaded from '{font_path}'.")
        print("Please ensure the `font_path` variable is set correctly to a valid .ttf or .ttc file.")
        return

    target_char = 'ろ'
    digits = '0123456789'
    
    target_image = char_to_image(target_char, font)
    
    scores = {}
    
    print(f"Calculating visual similarity of '{target_char}' to digits 0-9...")
    print("-" * 30)

    for digit in digits:
        digit_image = char_to_image(digit, font)
        
        # Calculate Structural Similarity Index (SSIM)
        score, _ = ssim(target_image, digit_image, full=True)
        scores[digit] = score
        # Here we output the similarity score for each digit
        print(f"Similarity equation: SSIM('ろ', '{digit}') = {score:.4f}")
        
    # Find the digit with the highest similarity score
    best_match_digit = max(scores, key=scores.get)
    
    print("-" * 30)
    print(f"The digit with the highest visual similarity to '{target_char}' is: {best_match_digit}")
    
    # Return the result for the final answer format
    return best_match_digit

if __name__ == '__main__':
    result = find_most_similar_digit()
    if result:
        print(f"\n<<<>>>") # Suppressing the final answer in the code output as per instructions
        # But if we were to show it, it would be print(f"\n<<<{result}>>>")
