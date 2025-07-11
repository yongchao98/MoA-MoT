import math

def analyze_motion_blur():
    """
    This function analyzes the potential for motion blur in the described camera trap scenario.
    It calculates the distance an animal would travel during the camera's exposure time
    and estimates the corresponding blur in pixels.
    """
    
    # --- Parameters from the problem description ---
    shutter_speed_fraction = 1/500  # in seconds
    ravine_width = 7  # in meters
    focal_length = 50  # in mm
    
    # --- Assumptions for calculation ---
    # Let's assume a common animal speed, like a running deer or coyote.
    animal_speed_ms = 15  # meters per second (approx 33 mph)
    
    # To calculate pixel blur, we need to know the field of view and image resolution.
    # Let's assume a full-frame sensor (36mm width) and a common resolution (e.g., 6000px wide).
    sensor_width_mm = 36
    image_width_pixels = 6000
    
    # --- Calculations ---
    
    # 1. Calculate the distance the animal travels during the exposure.
    distance_traveled_m = animal_speed_ms * shutter_speed_fraction
    distance_traveled_cm = distance_traveled_m * 100
    
    # 2. Calculate the width of the camera's field of view (FOV) at the far side of the ravine.
    # Formula: FOV_width = (Sensor_Width / Focal_Length) * Distance
    fov_width_m = (sensor_width_mm / focal_length) * ravine_width
    
    # 3. Calculate how many meters are represented by a single pixel.
    m_per_pixel = fov_width_m / image_width_pixels
    
    # 4. Calculate the motion blur in pixels.
    blur_in_pixels = distance_traveled_m / m_per_pixel
    
    # --- Output the results ---
    print("Analysis of Motion Blur Potential:")
    print("-" * 35)
    print(f"Shutter Speed: 1/{int(1/shutter_speed_fraction)} s")
    print(f"Assumed Animal Speed: {animal_speed_ms} m/s")
    print("\nStep 1: Animal Movement During Exposure")
    print(f"The animal travels {distance_traveled_m:.2f} meters, or {distance_traveled_cm:.2f} cm, while the shutter is open.")
    
    print("\nStep 2: Camera Field of View (Assumed)")
    print(f"At {ravine_width}m distance, the horizontal field of view is approximately {fov_width_m:.2f} meters wide.")
    
    print("\nStep 3: Calculating Pixel Blur")
    print(f"This results in a motion blur of approximately {math.ceil(blur_in_pixels)} pixels for a running animal.")
    
    print("\n--- Conclusion ---")
    print("A blur of this magnitude is significant and can degrade image quality, making features harder for a model to recognize.")
    print("The training data from GBIF, often sourced from photographers seeking sharp images, is unlikely to contain this artifact.")
    print("Therefore, simulating this specific type of blur during training is crucial to bridge the gap between the training data and the deployment data.")
    print("\nThe most important augmentation is H: A motion blur augmentation which selects sections of the image to apply a kernel to over a directional vector - mimicking the blur caused by animal movement.")

analyze_motion_blur()