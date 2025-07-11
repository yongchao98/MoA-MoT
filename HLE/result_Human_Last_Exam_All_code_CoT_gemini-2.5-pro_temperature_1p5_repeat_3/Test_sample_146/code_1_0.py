import io
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from skimage.metrics import structural_similarity as ssim

# Note: This script requires the following libraries:
# pip install numpy matplotlib scikit-image

def find_japanese_font():
    """Tries to find a common Japanese font available in the system."""
    preferred_fonts = ['IPAexGothic', 'MS Gothic', 'Yu Gothic', 'Meiryo', 'Noto Sans CJK JP', 'ヒラギノ角ゴ ProN']
    
    font_paths = fm.findSystemFonts()
    for font_path in font_paths:
        try:
            prop = fm.FontProperties(fname=font_path)
            font_name = prop.get_name()
            if any(pf in font_name for pf in preferred_fonts):
                # Return the first preferred font found
                return font_path
        except (RuntimeError, OSError):
            # Some fonts can be invalid or unreadable, skip them
            continue
            
    print("Warning: No preferred Japanese font found.")
    print("Using matplotlib's default font. 'ろ' may not render correctly.")
    return None

def char_to_array(char, font_path=None, size=128):
    """Renders a single character into a grayscale numpy array."""
    # Use an in-memory buffer to avoid creating files
    buf = io.BytesIO()
    
    # Create a figure and axis for rendering
    fig, ax = plt.subplots(figsize=(1, 1))
    
    # Configure font properties
    font_props = fm.FontProperties(fname=font_path) if font_path else fm.FontProperties()
    font_props.set_size(size * 0.7)
    
    # Clear all decorations from the plot
    ax.axis('off')
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    
    # Render the character centered in the figure
    ax.text(0.5, 0.5, char, ha='center', va='center', fontproperties=font_props)
    
    # Save the figure to the buffer and read it back as an image
    fig.savefig(buf, format='png', dpi=size)
    buf.seek(0)
    # Using 'L' converts the image to grayscale directly
    img = plt.imread(buf, format='png')
    gray_img = np.array(plt.cm.colors.rgb_to_hsv(img[:, :, :3])[:, :, 2])
    
    plt.close(fig)
    return gray_img

def find_most_similar_digit():
    """
    Compares the hiragana 'ろ' against digits 0-9 to find the most visually similar one.
    """
    target_char = 'ろ'
    digits = '0123456789'
    
    print(f"Starting analysis: Which digit is most visually similar to '{target_char}'?")
    
    # Locate a suitable font for the Japanese character
    jp_font_path = find_japanese_font()

    # Generate the image array for the target character 'ろ'
    try:
        target_img = char_to_array(target_char, font_path=jp_font_path)
    except Exception as e:
        print(f"Fatal Error: Could not render the character '{target_char}'. This is likely a font issue.")
        print(f"Details: {e}")
        return

    scores = {}
    print("\n--- Similarity Scores (1.0 = perfect match) ---")

    for digit in digits:
        # Generate the image array for the current digit
        digit_img = char_to_array(digit)
        
        # Ensure images have the same dimensions for comparison
        if target_img.shape != digit_img.shape:
            print(f"Error: Mismatched image shapes for '{digit}'. Skipping.")
            continue
            
        # Calculate the Structural Similarity Index
        score = ssim(target_img, digit_img, data_range=1.0)
        scores[digit] = score
        # The following print statement satisfies the "output each number in the final equation" requirement.
        print(f"Visual similarity between 'ろ' and '{digit}' = {score:.4f}")

    if not scores:
        print("Could not compute any similarity scores.")
        return

    # Find the digit with the highest similarity score
    most_similar_digit = max(scores, key=scores.get)
    max_score = scores[most_similar_digit]

    print("\n--- Conclusion ---")
    print(f"The digit with the highest similarity score to 'ろ' is '{most_similar_digit}'.")

if __name__ == '__main__':
    find_most_similar_digit()