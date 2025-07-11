import cv2
import numpy as np
import requests
from PIL import Image
import io

def solve_image_processing_task():
    """
    This script downloads an original and processed image, applies five different
    image processing techniques to the original, and calculates the Mean Squared Error (MSE)
    between each result and the target processed image to identify the correct technique.
    """

    # Function to download an image from a URL and return it in OpenCV BGR format
    def download_image(url):
        try:
            response = requests.get(url, timeout=15)
            response.raise_for_status()
            img_bytes = io.BytesIO(response.content)
            img_pil = Image.open(img_bytes).convert('RGB')
            return cv2.cvtColor(np.array(img_pil), cv2.COLOR_RGB2BGR)
        except requests.exceptions.RequestException as e:
            print(f"Error: Could not download image from {url}. {e}")
            return None

    # URLs for the original and target processed images
    # Using specific crops provided in the prompt for analysis
    url_original = "https://i.ibb.co/zXpXy1J/i-0-0.png"
    url_processed_target = "https://i.ibb.co/LgLdM3x/i-0-2.png"

    # Download images
    print("Downloading images...")
    original_img = download_image(url_original)
    target_img = download_image(url_processed_target)

    if original_img is None or target_img is None:
        print("Aborting analysis due to image download failure.")
        return

    print("Images downloaded successfully. Starting analysis...\n")
    h, w, _ = original_img.shape
    errors = {}

    def mean_squared_error(img1, img2):
        # Calculate MSE. The 'equation' is the mean of squared differences.
        # We will output the final result of this equation.
        return np.mean((img1.astype("float") - img2.astype("float")) ** 2)

    # --- Option A: Averaging filter, downsample, NN upscale ---
    img_a_blur = cv2.blur(original_img, (4, 4))
    img_a_down = cv2.resize(img_a_blur, (w // 4, h // 4), interpolation=cv2.INTER_NEAREST)
    img_a_final = cv2.resize(img_a_down, (w, h), interpolation=cv2.INTER_NEAREST)
    errors['A'] = mean_squared_error(target_img, img_a_final)

    # --- Option B: DCT, zero high-frequency, inverse DCT ---
    img_b_final = np.zeros_like(original_img, dtype=np.float32)
    for i in range(3):
        channel = original_img[:, :, i].astype(np.float32)
        dct_channel = cv2.dct(channel)
        rows, cols = dct_channel.shape
        # Zero out bottom-right half of coefficients
        dct_channel[int(rows * 0.75):, :] = 0
        dct_channel[:, int(cols * 0.75):] = 0
        img_b_final[:, :, i] = cv2.idct(dct_channel)
    img_b_final = np.clip(img_b_final, 0, 255).astype(np.uint8)
    errors['B'] = mean_squared_error(target_img, img_b_final)

    # --- Option C: Gaussian filter ---
    img_c_final = cv2.GaussianBlur(original_img, (7, 7), 0)
    errors['C'] = mean_squared_error(target_img, img_c_final)

    # --- Option D: Non-Local Means filter ---
    # Parameters from prompt: 7x7 template, 21x21 search window.
    # The filter strength 'h' is tuned to find the best match. A value around 10 is typical.
    img_d_final = cv2.fastNlMeansDenoisingColored(original_img, None, h=10, hColor=10, templateWindowSize=7, searchWindowSize=21)
    errors['D'] = mean_squared_error(target_img, img_d_final)

    # --- Option E: Downsample, Bilinear upscale ---
    img_e_down = cv2.resize(original_img, (w // 4, h // 4), interpolation=cv2.INTER_LINEAR)
    img_e_final = cv2.resize(img_e_down, (w, h), interpolation=cv2.INTER_LINEAR)
    errors['E'] = mean_squared_error(target_img, img_e_final)

    # --- Print results ---
    print("Results of Mean Squared Error (MSE) Calculation:")
    print("------------------------------------------------")
    print("The MSE measures the difference between the target image and the result of each option.")
    print("A lower MSE indicates a better match.\n")

    for option, error in sorted(errors.items(), key=lambda item: item[1]):
        print(f"Option {option}: MSE = {error:.2f}")

    best_option = min(errors, key=errors.get)
    print(f"\nConclusion: Option '{best_option}' has the lowest MSE, making it the most likely method used.")

# Run the analysis
solve_image_processing_task()