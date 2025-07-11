import requests
from PIL import Image
import io

def find_secret_word():
    """
    Finds a secret word encrypted in an image using LSB steganography.
    
    The plan is based on visual clues in the image:
    1.  An anime character points to the left border, indicating where the secret is hidden.
    2.  The left border is 8 pixels wide, hinting at 8-bit ASCII encoding.
    3.  The "noise" or "static" in the border is a classic sign of Least Significant Bit (LSB) steganography.
    
    The script will:
    - Fetch the image from a URL.
    - Read the pixels in the 8-pixel-wide left border.
    - Extract the least significant bit from the red channel of each pixel.
    - Group the bits into bytes and convert them to characters to reveal the word.
    """
    try:
        # The URL of the source image
        image_url = "https://i.imgur.com/vH3s11e.png"
        response = requests.get(image_url)
        # Raise an error if the image can't be downloaded
        response.raise_for_status()
        img_data = io.BytesIO(response.content)
        image = Image.open(img_data).convert('RGB')
    except Exception as e:
        print(f"Error: Could not retrieve or open the image. {e}")
        return

    pixels = image.load()
    width, height = image.size

    bits_string = ""
    border_width = 8
    # The top border is ~16 pixels high, so we start reading below it
    start_y = 16

    # Iterate through each pixel in the left border
    for y in range(start_y, height):
        for x in range(border_width):
            # Get the RGB value of the pixel
            r, g, b = pixels[x, y]
            # Extract the least significant bit (LSB) of the red channel
            # by performing a bitwise AND with 1.
            bits_string += str(r & 1)
            
    secret_message = ""
    # Process the bit string in 8-bit chunks (bytes)
    for i in range(0, len(bits_string), 8):
        byte = bits_string[i:i+8]
        if len(byte) == 8:
            # Convert the binary byte to an integer, then to an ASCII character
            char_code = int(byte, 2)
            secret_message += chr(char_code)
            
    # The full decoded message is often padded with null characters ('\x00')
    # at the end. We strip these to get the clean final word.
    final_word = secret_message.strip('\x00')
    
    print(f"The secret word hidden in the image is: {final_word}")

find_secret_word()