import cv2
import numpy as np

# --- User setup ---
# Please save the original image from the problem as 'original_image.png'
# in the same directory where you run this script.
input_filename = 'original_image.png'
output_filename = 'processed_with_nlm.png'

try:
    # Read the original image
    original_image = cv2.imread(input_filename)
    if original_image is None:
        raise FileNotFoundError

    # D. Apply a Non-Local Means filter with 7x7 pixel template patch and 21x21 pixel search window.
    # The 'h' parameter controls the filter strength. A higher 'h' means more smoothing.
    # We choose a relatively high value to match the strong smoothing in the example.
    h_parameter = 30 
    template_window_size = 7
    search_window_size = 21

    print(f"Applying Non-Local Means filter...")
    print(f"Template Patch Size: {template_window_size}x{template_window_size}")
    print(f"Search Window Size: {search_window_size}x{search_window_size}")
    
    # Apply the Non-Local Means filter for colored images
    processed_image = cv2.fastNlMeansDenoisingColored(
        original_image, 
        None, 
        h=h_parameter, 
        hColor=h_parameter, 
        templateWindowSize=template_window_size, 
        searchWindowSize=search_window_size
    )

    # Save the resulting image
    cv2.imwrite(output_filename, processed_image)
    print(f"\nProcessing complete.")
    print(f"The result has been saved as '{output_filename}'.")
    print("You can now visually compare this file with the 'Processed Image' from the problem.")

except FileNotFoundError:
    print(f"Error: The file '{input_filename}' was not found.")
    print("Please make sure you have saved the original image correctly in the same folder as the script.")
except Exception as e:
    print(f"An error occurred: {e}")
