def solve_talbot_history_problem():
    """
    This function provides the solution to the historical questions about the Talbot Settlement.
    The data is based on historical records.
    """
    
    # Historical data points
    # By 1823, the Talbot Settlement received a large influx of destitute migrants, particularly from 1820-1823.
    destitute_migrants = 12000
    
    # Colonel Talbot's initial personal land grant in 1803.
    initial_grant_acres = 5000
    
    # The total acreage he eventually came to administer and claim authority over.
    total_claimed_acres = 650000
    
    # Calculate how much larger the claimed acreage was than the original grant.
    acreage_difference = total_claimed_acres - initial_grant_acres
    
    # Print the answers
    print(f"As a result of the land grant, approximately {destitute_migrants} destitute migrants settled during the peak years of 1820-1823.")
    print("The acreage he eventually claimed was larger than the original land grant by this many acres:")
    print(f"{total_claimed_acres} (total claimed) - {initial_grant_acres} (initial grant) = {acreage_difference} acres")

solve_talbot_history_problem()