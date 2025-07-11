import sys
from PIL import Image
import os
import re

def find_secret():
    """
    This function analyzes the provided image to find a secret message,
    which in turn reveals a secret word and an equation.
    The message is encoded in the pixel colors of the left border.
    """
    
    # The image is expected to be in the same directory as the script.
    image_path = 'image.png'
    
    if not os.path.exists(image_path):
        print(f"Error: The image file '{image_path}' was not found.")
        print("Please ensure the image is in the correct directory.")
        return

    try:
        img = Image.open(image_path).convert('RGB')
        pixels = img.load()
    except Exception as e:
        print(f"Error opening or processing image: {e}")
        return

    def get_bit_triplet_from_color(rgb_tuple):
        """
        Converts an RGB color tuple into a 3-bit binary string.
        A color channel value > 128 is considered a '1', otherwise '0'.
        The order is R, G, B for the bits from most to least significant.
        """
        r, g, b = rgb_tuple
        bit2 = '1' if r > 128 else '0'
        bit1 = '1' if g > 128 else '0'
        bit0 = '1' if b > 128 else '0'
        return f"{bit2}{bit1}{bit0}"

    # The data is in an 8-pixel wide vertical strip on the left.
    # We sample from the middle of this strip (x=12).
    # The data blocks are 8 pixels high and start at y=40.
    # We sample from the middle of each block (y = 44, 52, 60, ...).
    x_coord = 12
    y_start = 44
    y_step = 8
    
    # The main content area where the border exists is 504 pixels high (from y=40 to 544).
    # This corresponds to 504 / 8 = 63 data blocks.
    num_blocks = 63
    
    bit_stream = ""
    for i in range(num_blocks):
        y = y_start + i * y_step
        pixel_color = pixels[x_coord, y]
        bit_triplet = get_bit_triplet_from_color(pixel_color)
        bit_stream += bit_triplet
        
    # The total bit stream has 63 * 3 = 189 bits.
    # We group these bits into 8-bit bytes to form ASCII characters.
    # 189 bits will form 23 full bytes with 5 bits remaining, which we'll ignore.
    byte_values = []
    for i in range(0, len(bit_stream) // 8):
        byte_str = bit_stream[i*8 : (i+1)*8]
        byte_values.append(int(byte_str, 2))

    # The puzzle uses the first two bytes as numbers for an equation.
    if len(byte_values) < 3:
        print("Not enough data decoded to find the equation.")
        return
        
    num1 = byte_values[0]
    num2 = byte_values[1]
    
    # The puzzle's solution is the equation derived from these numbers.
    result = num1 - num2
    
    # Print the equation, showing each number as requested by the prompt.
    print(f"{num1} - {num2} = {result}")

# Run the solver function
find_secret()