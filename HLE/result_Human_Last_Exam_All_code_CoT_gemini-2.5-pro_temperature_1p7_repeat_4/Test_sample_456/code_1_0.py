import cv2
import albumentations as A
import numpy as np

def run_augmentation_demo():
    """
    This function demonstrates the effect of motion blur augmentation.
    It creates a sample image with a shape and applies motion blur to it.
    This simulates how an augmentation pipeline would prepare training data
    to make a model robust to motion artifacts expected in the deployment scenario.
    """

    # 1. Create a simple dummy image (e.g., 256x256 with a white circle)
    # This represents a clear photo of an animal from the training set.
    image = np.zeros((256, 256, 3), dtype=np.uint8)
    cv2.circle(image, (128, 128), 50, (255, 255, 255), -1)

    # 2. Define the motion blur augmentation.
    # The 'blur_limit' parameter controls the severity of the blur.
    # The 'p=1.0' means the augmentation will always be applied for this demo.
    transform = A.Compose([
        A.MotionBlur(blur_limit=(21, 41), p=1.0),
    ])

    # 3. Apply the augmentation to the image
    transformed = transform(image=image)
    blurred_image = transformed['image']

    # 4. Print explanation and show results (in a real scenario, we wouldn't show)
    # In a real environment, you might just get a file path or nothing, 
    # but here we print a confirmation.
    print("Original Image: A clear, static shape.")
    print("Augmented Image: The shape is blurred horizontally, mimicking movement.")
    print("\nExplanation:")
    print("The most significant difference between the curated training data (e.g., from GBIF) and the camera trap's deployment data is that the deployment images will capture animals in motion at night.")
    print("This movement creates motion blur, which significantly changes the animal's apparent shape and texture.")
    print("Training the model with a Motion Blur augmentation (H) teaches it to recognize species even when their features are distorted by movement.")
    print("This is the most critical augmentation for ensuring accuracy in the specific deployment environment described.")

# Execute the demonstration
if __name__ == '__main__':
    run_augmentation_demo()
