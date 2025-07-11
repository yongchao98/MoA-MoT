import requests
import io
import math
from PIL import Image

def solve():
    """
    Downloads the image, extracts the hidden message from the left border,
    decodes it, and prints the secret word.
    """
    # The URL of the image containing the secret message
    image_url = "https://storage.googleapis.com/testmaker-prod-project-184518.appspot.com/assets/images/30372_1_400_600_400_600.png"

    # Step 1: Define the mapping from colors to 3-bit values (trits)
    # The value is calculated as: 4*Blue_bit + 2*Green_bit + 1*Red_bit
    COLOR_MAP = {
        (0, 0, 0): 0,        # Black: 000
        (255, 0, 0): 1,      # Red:   001
        (0, 255, 0): 2,      # Green: 010
        (255, 255, 0): 3,    # Yellow:011
        (0, 0, 255): 4,      # Blue:  100
        (255, 0, 255): 5,    # Magenta:101
        (0, 255, 255): 6,    # Cyan:  110
        (255, 255, 255): 7   # White: 111
    }

    # Helper function to find the closest color in our map for a given pixel
    def get_closest_color_value(pixel_rgb):
        min_dist = float('inf')
        value = -1
        # Unpack the pixel's RGB, ignoring alpha if it exists
        r1, g1, b1 = pixel_rgb[:3]
        for color_rgb, trit_val in COLOR_MAP.items():
            r2, g2, b2 = color_rgb
            # Calculate Euclidean distance in RGB space
            dist = math.sqrt((r1 - r2)**2 + (g1 - g2)**2 + (b1 - b2)**2)
            if dist < min_dist:
                min_dist = dist
                value = trit_val
        return value

    # Step 2: Download and load the image
    try:
        response = requests.get(image_url)
        response.raise_for_status()  # Raise an exception for bad status codes
        img = Image.open(io.BytesIO(response.content))
        pixels = img.load()
    except Exception as e:
        print(f"Error loading image: {e}")
        return

    # Step 3: Extract the bitstream from the left border
    bit_stream = []
    # The data is in a column on the left border. We'll use x=3.
    # The top and bottom borders are 16px high, so we scan from y=16 to y=583.
    x_coord = 3
    for y in range(16, 584):
        pixel_rgb = pixels[x_coord, y]
        trit_value = get_closest_color_value(pixel_rgb)
        
        # Convert the 3-bit value (0-7) to three individual bits and add to the stream
        # (MSB first, corresponding to Blue, Green, Red bits)
        bit_stream.append((trit_value >> 2) & 1)
        bit_stream.append((trit_value >> 1) & 1)
        bit_stream.append(trit_value & 1)

    # Step 4: Convert the bitstream to ASCII text
    secret_message = ""
    # Group the bitstream into 8-bit chunks (bytes)
    for i in range(0, len(bit_stream), 8):
        byte_bits = bit_stream[i:i+8]
        if len(byte_bits) == 8:
            # Join bits into a string '10101010', convert to integer, then to character
            byte_str = "".join(map(str, byte_bits))
            char_code = int(byte_str, 2)
            secret_message += chr(char_code)

    # Step 5: Print the final decoded message
    print(secret_message)

solve()