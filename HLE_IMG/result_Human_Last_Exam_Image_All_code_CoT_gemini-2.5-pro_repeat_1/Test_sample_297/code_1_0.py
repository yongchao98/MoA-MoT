import requests
from PIL import Image
from io import BytesIO
import math

def solve_secret_word():
    """
    Finds the secret word encoded in the image's right border.
    The method involves:
    1. Fetching the image from a URL.
    2. Sampling the 100 color blocks from the right border.
    3. Mapping each of the 8 possible colors to a 3-bit value using a BGR (Blue, Green, Red) scheme.
    4. Concatenating the bits to form a single bitstream.
    5. Grouping the bitstream into 7-bit chunks.
    6. Converting each chunk from a binary string to an ASCII character.
    7. Joining the characters to reveal the secret word.
    """
    try:
        # Step 1: Fetch the image
        url = "https://i.imgur.com/u1GUOcr.jpg"
        response = requests.get(url)
        response.raise_for_status()
        img = Image.open(BytesIO(response.content)).convert('RGB')
    except requests.exceptions.RequestException as e:
        print(f"Error fetching image: {e}")
        return
    except Exception as e:
        print(f"Error processing image: {e}")
        return

    width, height = img.size

    # The 8 pure colors used in the border
    PURE_COLORS = [
        (0, 0, 0), (255, 0, 0), (0, 255, 0), (0, 0, 255),
        (255, 255, 0), (255, 0, 255), (0, 255, 255), (255, 255, 255)
    ]

    def find_closest_color(rgb):
        """Finds the closest pure color to a given RGB tuple."""
        min_dist = float('inf')
        closest_color = None
        for color in PURE_COLORS:
            dist = math.sqrt(sum([(c1 - c2) ** 2 for c1, c2 in zip(rgb, color)]))
            if dist < min_dist:
                min_dist = dist
                closest_color = color
        return closest_color

    # Step 2: Sample the color blocks
    color_sequence = []
    x_coord = width - 6  # A column within the right border
    # The color blocks are 8px high and start at y=28.
    # We sample the center of each block.
    for y in range(32, 829, 8):
        pixel_rgb = img.getpixel((x_coord, y))
        closest_pure_color = find_closest_color(pixel_rgb)
        color_sequence.append(closest_pure_color)
    
    # Step 3: Define the BGR color-to-bit mapping
    # val = Blue*4 + Green*2 + Red*1
    def color_to_bits(color):
        r = 1 if color[0] > 128 else 0
        g = 1 if color[1] > 128 else 0
        b = 1 if color[2] > 128 else 0
        value = b * 4 + g * 2 + r * 1
        return format(value, '03b')

    # Step 4: Create the bitstream
    bit_stream = "".join([color_to_bits(c) for c in color_sequence])

    # Step 5 & 6: Group by 7 bits and convert to ASCII
    secret_message = []
    for i in range(0, len(bit_stream), 7):
        chunk = bit_stream[i:i+7]
        if len(chunk) == 7:
            char_code = int(chunk, 2)
            # Only consider printable ASCII characters
            if 32 <= char_code <= 126:
                secret_message.append(chr(char_code))

    # Step 7: Print the final result
    result = "".join(secret_message)
    print(f"The secret word is: {result}")

solve_secret_word()