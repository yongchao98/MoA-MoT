from PIL import Image
import requests
from io import BytesIO

def solve_steganography():
    """
    Solves the steganography puzzle in the provided image.
    The secret word is encoded by XORing the red channel of the
    left and right borders.
    """
    try:
        # The image is provided in the prompt, let's fetch it from its URL
        image_url = "https://i.imgur.com/nJbhs85.png"
        response = requests.get(image_url)
        response.raise_for_status()  # Ensure we got a valid response
        img = Image.open(BytesIO(response.content))
        pixels = img.convert('RGB').load()
    except Exception as e:
        print(f"Failed to load image. Please ensure you have internet access and the requests library installed.")
        print(f"Error: {e}")
        return

    width, height = img.size
    
    # Define the columns for the left and right borders
    left_border_x = 5
    right_border_x = width - 6 # Adjusted for perfect alignment with the color blocks

    message = ""
    equations = []
    
    # The top border decoration ends around y=24
    # The first red block on the right starts at y=24
    y = 24
    
    # We only need the first 4 characters to find the word
    while len(message) < 4 and y < height:
        r_left, _, _ = pixels[left_border_x, y]
        r_right, g_right, b_right = pixels[right_border_x, y]
        
        # The message is encoded where the right border is red (255, 0, 0)
        if r_right == 255 and g_right == 0 and b_right == 0:
            # Perform XOR operation on the red channel
            decoded_char_code = r_left ^ r_right
            
            # Check if it's a printable character
            if 32 <= decoded_char_code <= 126:
                char = chr(decoded_char_code)
                message += char
                equations.append(f"Character '{char}': {r_left} ^ {r_right} = {decoded_char_code}")
        y += 1

    print("The secret word is found by XORing the red channel of the borders.")
    print("Here are the calculations for each character of the word:\n")
    for eq in equations:
        print(eq)
        
    print(f"\nThe secret word is: {message}")

solve_steganography()