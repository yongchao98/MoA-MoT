import requests
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from skimage.restoration import denoise_nl_means, estimate_sigma
from skimage.util import img_as_float
from io import BytesIO

def apply_non_local_means_filter():
    """
    Downloads the parrot image, applies a Non-Local Means filter,
    and displays the original and processed images side-by-side.
    """
    try:
        # URL of the original image
        url = "https://i.imgur.com/z4b0yv2.png"
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        # Open the image and convert to a numpy array of floats
        original_img = Image.open(BytesIO(response.content))
        # The provided image is PNG and might have an alpha channel, remove it.
        if original_img.mode == 'RGBA':
            original_img = original_img.convert('RGB')
            
        img = img_as_float(original_img)

        # Parameters for the NLM filter as described in Option D
        # A 7x7 pixel template patch
        patch_size = 7
        # A 21x21 pixel search window. patch_distance is the radius of the search window.
        # So, for a 21x21 window, the radius is (21-1)/2 = 10.
        patch_distance = 10
        
        # The h parameter controls filtering strength. It's often related to noise sigma.
        # We will set a value that produces a visually strong smoothing effect.
        # A larger h leads to more smoothing.
        h_value = 0.08
        
        print(f"Applying Non-Local Means filter...")
        print(f"Parameters: patch_size={patch_size}, patch_distance={patch_distance}, h={h_value}")

        # Apply the Non-Local Means filter
        # channel_axis=-1 is used for multichannel (RGB) images
        denoised_img = denoise_nl_means(
            img,
            h=h_value,
            patch_size=patch_size,
            patch_distance=patch_distance,
            channel_axis=-1,
            fast_mode=True
        )
        
        # Display the results
        fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharex=True, sharey=True)
        ax = axes.ravel()

        ax[0].imshow(img)
        ax[0].set_title('Original Image')
        ax[0].axis('off')

        ax[1].imshow(denoised_img)
        ax[1].set_title('Processed with Non-Local Means')
        ax[1].axis('off')

        plt.tight_layout()
        plt.show()

    except requests.exceptions.RequestException as e:
        print(f"Error downloading the image: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have scikit-image, numpy, matplotlib, and requests installed.")
        print("You can install them using: pip install scikit-image numpy matplotlib requests Pillow")

# Run the function
apply_non_local_means_filter()
