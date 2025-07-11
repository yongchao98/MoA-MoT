def describe_azobenzene_table():
    """
    This function describes the functionally accurate behavior of a picnic table
    shaped like an azobenzene molecule throughout a day.
    """
    
    sunrise_description = (
        "At sunrise, when exposed to the sun's UV radiation, the table "
        "must change its shape. It should convert from its stable, elongated 'trans' form "
        "into its less stable, bent 'cis' form."
    )
    
    sunset_description = (
        "At sunset, as the UV light source disappears, the table must "
        "revert to its original state. The bent 'cis' form should relax back "
        "into the stable, elongated 'trans' form, ready for the next day's cycle."
    )
    
    print("To be functionally accurate, the azobenzene-shaped picnic table must do the following:")
    print("\n[ACTION AT SUNRISE]")
    print(sunrise_description)
    print("\n[ACTION AT SUNSET]")
    print(sunset_description)

if __name__ == '__main__':
    describe_azobenzene_table()