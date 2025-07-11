def generate_surveillance_plan():
    """
    This function outlines the appropriate post-procedure surveillance program
    based on clinical best practices for SFA stenting.
    """
    
    # Define the key components of the correct surveillance plan
    base_assessment = "Regular visits with assessment for interval change in symptoms, vascular examination"
    imaging_modality = "arterial duplex"
    time_points = [3, 6, 12]
    long_term_follow_up_years = 2

    # Construct the descriptive sentence for the plan
    # This ensures all numbers from the chosen option are part of the output
    plan_description = (
        f"{base_assessment}, and {imaging_modality} at "
        f"{time_points[0]} months, "
        f"{time_points[1]} months, "
        f"{time_points[2]} months, and "
        f"{long_term_follow_up_years} years"
    )

    print("The most appropriate post-procedure surveillance program is:")
    print(plan_description)

# Execute the function to print the answer
generate_surveillance_plan()