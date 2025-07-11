import cv2
import numpy as np
import matplotlib.pyplot as plt
import requests
from PIL import Image
from io import BytesIO

# --- Image Loading ---
# The user should provide the image file.
# As I don't have access to the local file system,
# I'll download the image from the provided URL.
# In a local environment, you would use:
# original_image = cv2.imread('parrot_original.png')
# original_image = cv2.cvtColor(original_image, cv2.COLOR_BGR2RGB)

try:
    # URL of the original image provided in the problem description
    image_url = "https://i.imgur.com/kSXYDZb.png"
    response = requests.get(image_url)
    response.raise_for_status()  # Raise an exception for bad status codes
    img = Image.open(BytesIO(response.content))
    original_image = np.array(img)
    # The image is loaded as RGBA, let's convert to RGB for processing
    if original_image.shape[2] == 4:
        original_image = original_image[:, :, :3]

    # --- Image Processing: Applying Non-Local Means Filter ---
    # Parameters from option D:
    # patch size = 7x7 pixels -> templateWindowSize = 7
    # search window = 21x21 pixels -> searchWindowSize = 21
    # 'h' is the filter strength parameter. A value of 10 is common.
    # 'hColor' is the strength for color components.
    processed_image = cv2.fastNlMeansDenoisingColored(
        original_image,
        None,
        h=10,
        hColor=10,
        templateWindowSize=7,
        searchWindowSize=21
    )

    # --- Displaying the result ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    axes[0].imshow(original_image)
    axes[0].set_title('Original Image')
    axes[0].axis('off')

    axes[1].imshow(processed_image)
    axes[1].set_title('Processed with Non-Local Means')
    axes[1].axis('off')

    plt.tight_layout()
    plt.show()

except requests.exceptions.RequestException as e:
    print(f"Failed to download the image. Please download it manually from 'https://i.imgur.com/kSXYDZb.png' and save it as 'parrot_original.png'.")
    print(f"Error: {e}")
except Exception as e:
    print(f"An error occurred: {e}")
    print("Please ensure the image file 'parrot_original.png' is in the same directory as the script and try again.")
