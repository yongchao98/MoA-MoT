def find_tornadogenesis_inhibitor():
    """
    Analyzes meteorological data to find the number with the greatest negative
    impact on tornadogenesis.
    """
    
    # The K-Index is a measure of thunderstorm potential. A value of 20 is very low
    # and suggests only isolated thunderstorm activity is likely. In an otherwise
    # explosive environment (indicated by high CAPE, high shear, and very high
    # STP/SRH values), this low potential for storm initiation is the most
    # significant limiting factor for tornado development. It indicates a potential
    # for the forecast to be a "bust" if robust storms cannot form.
    inhibiting_factor_value = 20
    inhibiting_factor_name = "K-Index"

    print(f"The number in the dataset that most significantly reduces the likelihood of tornado formation is from the {inhibiting_factor_name}.")
    print("This value indicates a very low potential for widespread thunderstorm development, which is a primary prerequisite for tornadogenesis.")
    print("Despite other parameters being extremely favorable, if storms cannot initiate, tornadoes cannot form.")
    print("\nInhibiting Factor:")
    print(f"{inhibiting_factor_name} = {inhibiting_factor_value}")

find_tornadogenesis_inhibitor()