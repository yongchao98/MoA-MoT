def surveillance_plan():
    """
    This function outlines the components of the most appropriate post-procedure
    surveillance program for SFA stenting based on clinical best practices.
    """
    
    # The optimal surveillance program uses the best imaging modality at a frequency
    # that matches the period of highest risk for restenosis.

    # 1. The preferred surveillance tool
    primary_tool = "Arterial Duplex Ultrasound"

    # 2. The surveillance schedule, with increased frequency in the first year.
    # The numbers represent the time points for follow-up visits with duplex ultrasound.
    first_year_follow_up_months = [3, 6, 12]
    second_year_follow_up_years = 2

    print("Recommended Post-Procedure Surveillance Program:")
    print("---------------------------------------------")
    print(f"Primary Surveillance Modality: {primary_tool}")
    print("\nSurveillance Schedule:")
    print(f"- At {first_year_follow_up_months[0]} months post-procedure")
    print(f"- At {first_year_follow_up_months[1]} months post-procedure")
    print(f"- At {first_year_follow_up_months[2]} months post-procedure")
    print(f"- At {second_year_follow_up_years} years post-procedure")
    print("\nThis plan, along with regular clinical assessment, is designed to detect in-stent restenosis early.")

surveillance_plan()