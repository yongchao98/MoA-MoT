import requests
from PIL import Image
from io import BytesIO

def find_secret_word():
    """
    Solves the puzzle by decoding a secret message from the image.
    The secret message is encrypted in the vertical color bar on the right side.
    
    The decoding process is as follows:
    1. The image is downloaded from its URL.
    2. Each color in the bar corresponds to a 3-bit binary value (e.g., Red -> '100').
    3. The script reads these colors from top to bottom, building a long binary string.
    4. This binary string is then chunked into 8-bit bytes.
    5. Each byte is converted to its corresponding ASCII character to reveal the secret word.
    """
    
    # URL of the image to be processed
    image_url = "https://i.imgur.com/epCgQ.png"
    
    try:
        # Step 1: Download the image from the URL
        response = requests.get(image_url)
        response.raise_for_status()  # Ensure the download was successful
        
        # Open the image using the Pillow library
        img = Image.open(BytesIO(response.content))
        
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the image. {e}")
        return
    except Exception as e:
        print(f"Error: Could not process the image. {e}")
        return

    # Step 2: Define the mapping from RGB colors to 3-bit binary strings
    COLOR_TO_BITS = {
        (255, 0, 0): '100',      # Red
        (0, 255, 0): '010',      # Green
        (0, 0, 255): '001',      # Blue
        (255, 255, 0): '110',    # Yellow
        (0, 255, 255): '011',    # Cyan
        (255, 0, 255): '101',    # Magenta
        (255, 255, 255): '111',  # White
        (0, 0, 0): '000',        # Black
    }

    pixels = img.convert('RGB').load()
    width, height = img.size
    
    # The color bar is consistently located 5 pixels from the right edge
    x_coord = width - 5
    binary_string = ""
    processed_y_coords = set()

    # Step 3: Iterate through the pixels of the color bar
    for y in range(height):
        # Skip this pixel if its color block has already been processed
        if y in processed_y_coords:
            continue
            
        pixel_color = pixels[x_coord, y]
        
        # Process only the key colors, ignoring black separators
        if pixel_color in COLOR_TO_BITS and pixel_color != (0, 0, 0):
            # Append the corresponding 3-bit string
            binary_string += COLOR_TO_BITS[pixel_color]
            
            # Mark all pixels in this colored block as processed to avoid duplicates
            y_block = y
            while y_block < height and pixels[x_coord, y_block] == pixel_color:
                processed_y_coords.add(y_block)
                y_block += 1
                
    # Step 4: Convert the binary string to ASCII text
    secret_message = ""
    for i in range(0, len(binary_string), 8):
        byte = binary_string[i:i+8]
        if len(byte) == 8:
            # Convert 8-bit string to an integer, then to a character
            decimal_value = int(byte, 2)
            secret_message += chr(decimal_value)
            
    # Step 5: Print the final result
    print(f"The decoded secret word is: {secret_message}")

# Execute the function to find the word
find_secret_word()