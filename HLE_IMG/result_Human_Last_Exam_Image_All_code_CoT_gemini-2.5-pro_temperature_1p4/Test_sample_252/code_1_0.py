import cv2
import numpy as np
import requests
from PIL import Image
import io

def analyze_image_processing():
    """
    Downloads the image, applies all 5 processing methods,
    and calculates the Mean Squared Error (MSE) for each
    to determine the most likely method.
    """
    try:
        # Step 1: Download and prepare images
        # The URL points to the collage of Original | Processed
        url = "https://i.imgur.com/G5g20iW.png"
        print(f"Downloading image from {url}...")
        response = requests.get(url, timeout=15)
        response.raise_for_status()  # Raise an exception for bad status codes
        collage_img_pil = Image.open(io.BytesIO(response.content))
        print("Image downloaded successfully.")

        # The collage is 1414x403. Left is original, right is processed.
        w, h = collage_img_pil.size
        mid_w = w // 2

        # Crop images from the collage
        original_pil = collage_img_pil.crop((0, 0, mid_w, h))
        processed_pil = collage_img_pil.crop((mid_w, 0, w, h))

        # Convert PIL images to OpenCV format (BGR) for processing
        original_image = cv2.cvtColor(np.array(original_pil), cv2.COLOR_RGB2BGR)
        ground_truth_processed = cv2.cvtColor(np.array(processed_pil), cv2.COLOR_RGB2BGR)
        
        img_h, img_w, _ = original_image.shape

    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the image. {e}")
        print("Please check the network connection and the URL.")
        return
    except Exception as e:
        print(f"An error occurred during image preparation: {e}")
        return

    # Step 2: Define MSE function
    def calculate_mse(imageA, imageB):
        # Calculate the Mean Squared Error between two images
        err = np.sum((imageA.astype("float") - imageB.astype("float")) ** 2)
        err /= float(imageA.shape[0] * imageA.shape[1] * imageA.shape[2])
        return err

    print("\nAnalyzing the image processing methods by comparing their output to the processed image...")
    print("---------------------------------------------------------")

    # --- Step 3: Apply and evaluate each option ---

    # Option A: Averaging filter (4x4), downsample by 4, upscale (nearest)
    avg_img_A = cv2.blur(original_image, (4, 4))
    down_A = cv2.resize(avg_img_A, (img_w // 4, img_h // 4), interpolation=cv2.INTER_AREA)
    res_A = cv2.resize(down_A, (img_w, img_h), interpolation=cv2.INTER_NEAREST)
    mse_A = calculate_mse(res_A, ground_truth_processed)
    print(f"Option A MSE (Avg Filter -> Downsample -> Upscale Nearest): {mse_A:.2f}")

    # Option B: DCT transform, zeroing half of the high frequencies
    res_B_float = original_image.copy().astype(np.float32)
    for i in range(3):
        channel = res_B_float[:, :, i]
        dct_channel = cv2.dct(channel)
        rows, cols = dct_channel.shape
        dct_channel[int(rows * 0.5):, :] = 0
        dct_channel[:, int(cols * 0.5):] = 0
        res_B_float[:, :, i] = cv2.idct(dct_channel)
    res_B = np.clip(res_B_float, 0, 255).astype(np.uint8)
    mse_B = calculate_mse(res_B, ground_truth_processed)
    print(f"Option B MSE (DCT High-Frequency Cutoff): {mse_B:.2f}")

    # Option C: Gaussian filter with a 7x7 kernel
    res_C = cv2.GaussianBlur(original_image, (7, 7), 0)
    mse_C = calculate_mse(res_C, ground_truth_processed)
    print(f"Option C MSE (Gaussian Filter 7x7): {mse_C:.2f}")

    # Option D: Non-Local Means filter with 7x7 patch and 21x21 window
    # h=filter strength, a value of 10 is standard for colored images.
    res_D = cv2.fastNlMeansDenoisingColored(original_image, None, h=10, hColor=10, templateWindowSize=7, searchWindowSize=21)
    mse_D = calculate_mse(res_D, ground_truth_processed)
    print(f"Option D MSE (Non-Local Means Filter): {mse_D:.2f}")

    # Option E: Downsample by 4, upscale with a bilinear filter
    down_E = cv2.resize(original_image, (img_w // 4, img_h // 4), interpolation=cv2.INTER_AREA)
    res_E = cv2.resize(down_E, (img_w, img_h), interpolation=cv2.INTER_LINEAR)
    mse_E = calculate_mse(res_E, ground_truth_processed)
    print(f"Option E MSE (Downsample -> Upscale Bilinear): {mse_E:.2f}")

    print("---------------------------------------------------------")

    # Step 4: Find and print the best option
    mse_values = {'A': mse_A, 'B': mse_B, 'C': mse_C, 'D': mse_D, 'E': mse_E}
    best_option = min(mse_values, key=mse_values.get)
    print(f"\nConclusion: Option '{best_option}' has the lowest Mean Squared Error.")
    print("This indicates it is the most likely processing method used.")

if __name__ == '__main__':
    analyze_image_processing()