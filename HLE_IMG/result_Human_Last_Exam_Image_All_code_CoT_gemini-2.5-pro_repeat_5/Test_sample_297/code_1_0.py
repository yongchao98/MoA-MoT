import requests
from PIL import Image
import io

def find_secret_word():
    """
    Downloads the image, extracts the color data from the right border,
    and decodes the hidden message.
    """
    try:
        # 1. Download and open the image
        url = "https://i.imgur.com/Gj3oPz3.png"
        response = requests.get(url)
        response.raise_for_status()
        img = Image.open(io.BytesIO(response.content)).convert('RGB')
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the image. {e}")
        return
    except Exception as e:
        print(f"Error: Could not process the image. {e}")
        return

    width, height = img.size

    # 2. Define the mapping from colors to their 3-bit binary representations
    color_to_bits = {
        (0, 0, 0): '000',     # Black
        (255, 0, 0): '001',   # Red
        (0, 255, 0): '010',   # Green
        (255, 255, 0): '011', # Yellow
        (0, 0, 255): '100',   # Blue
        (255, 0, 255): '101', # Magenta
        (0, 255, 255): '110', # Cyan
        (255, 255, 255): '111'  # White
    }

    # 3. Extract the sequence of colors from the vertical bar on the right
    x_coord = width - 5  # Sample from a column safely within the border
    extracted_colors = []
    last_color_found = None
    
    for y in range(height):
        pixel_rgb = img.getpixel((x_coord, y))
        
        # Find the closest color in our palette to handle compression artifacts
        min_dist = float('inf')
        closest_palette_color = None
        for p_rgb in color_to_bits.keys():
            dist = sum((a - b)**2 for a, b in zip(pixel_rgb, p_rgb))
            if dist < min_dist:
                min_dist = dist
                closest_palette_color = p_rgb

        # Add the color to our list if it's a new, distinct block
        if closest_palette_color != last_color_found and min_dist < 5000:
             # Skip the pure black background at the very top of the image
            if y > 10:
                extracted_colors.append(closest_palette_color)
            last_color_found = closest_palette_color

    # 4. Identify the message data, which comes after a header sequence
    header_sequence = [
        (255, 0, 0), (0, 255, 0), (0, 0, 255), (255, 255, 0),
        (0, 255, 255), (255, 0, 255), (255, 255, 255)
    ]
    
    message_start_index = 0
    # The header might appear multiple times; we want the data after the last one
    for i in range(len(extracted_colors) - len(header_sequence)):
        if extracted_colors[i:i+len(header_sequence)] == header_sequence:
            message_start_index = i + len(header_sequence)

    message_colors = extracted_colors[message_start_index:]

    # 5. Convert the sequence of message colors into a single bit stream
    bit_stream = "".join([color_to_bits[c] for c in message_colors])

    # 6. Decode the bit stream into ASCII characters
    secret_word = ""
    # Process the stream in 8-bit chunks (bytes)
    for i in range(0, len(bit_stream) - len(bit_stream) % 8, 8):
        byte = bit_stream[i:i+8]
        char_code = int(byte, 2)
        secret_word += chr(char_code)

    print(f"The decoded secret word is: {secret_word}")

if __name__ == "__main__":
    find_secret_word()